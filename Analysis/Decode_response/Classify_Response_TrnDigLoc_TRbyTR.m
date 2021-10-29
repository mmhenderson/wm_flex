%% Response decoding analysis
% Train and test linear decoder, using data from one task condition at a
% time. Labels are the expected (correct) response on each trial. Training
% and testing on each timepoint for time-resolved decoding (TR by TR)
% cross-validate across sessions, response/luminance mapping was swapped
% between sessions.
% Saves the results in a mat file, can plot it using a separate script
% (plotClassResults_Response_TRbyTR.m)
%%
clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;
filesepinds = find(curr_dir==filesep);
nDirsUp = 2;
exp_path = curr_dir(1:filesepinds(end-nDirsUp+1));
addpath(fullfile(exp_path,'Analysis','stats_code'));

nVox2Use = 10000;    % this is the max number of vox to use, so if it's very big we're using all the voxels.
nPermIter = 1000;       % for generating null decoding accuracies, how many iterations of shuffling to do?

% what kind of classifier using?
class_str = 'normEucDist';
% get ready for parallel pool operations
dbstop if error
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 534544;
rng(rndseed,'twister');

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);
nTrialsTotal=2*10*20;

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    fn2load = fullfile(exp_path,'Samples',sprintf('DigLocSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponse_TrnDigLoc_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

    v2do=[1:length(mainSig)];
    areas2test = [1:length(mainSig)];

    for vv = v2do
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        %% loop over conditions
        for cc =1:nConds

            %% first organize my training data (digit localizer)
            
            % using all trials because accuracy is generally very high. Can
            % also remove the incorrect trials, but that would un-balance
            % the training set.
            trials2use = ones(size(locSig(1).ExpDigit,1),1)==1;
%             trials2use = locSig(vv).ActDigit == locSig(vv).ExpDigit;
%             trnLabs = locSig(vv).ActDigit(trials2use);
            % using "expected" digit here, to keep training set balanced.
            trnLabs = locSig(vv).ExpDigit(trials2use);
            assert(sum(trnLabs==1)==sum(trnLabs==2));
            
%             trnRuns = locSig(vv).runLabs(trials2use);
            trnDat = locSig(vv).dat_avg(trials2use,:);
            trnDat = trnDat - repmat(mean(trnDat,2), 1, size(trnDat,2));
            
            %% now organize my testing data (main task)
            
            tstLabs = mainSig(vv).CorrectResp;
            condLabsTest = mainSig(vv).condLabs;
%             tstRuns = mainSig(vv).runLabs;
            % getting rid of any trials with no response here
            % also taking out just the relevant condition!!
            trials2use = condLabsTest==cc;
            
            tstLabs = tstLabs(trials2use);

            % dat is [ntrials x nTRs x nVox]
            mainDat = mainSig(vv).dat_by_TR;
            mainDat = mainDat(trials2use,:,:);

            % subtract mean over voxels 
            mainDat = mainDat - repmat(mean(mainDat,3), 1, 1, size(mainDat, 3));

            % mainDat is [ntrials x nTRs x nVox]
            nTRs_out = size(mainDat,2);
 
            if vv==v2do(1) && cc==1
                % preallocate array here
                allacc = nan(length(ROI_names), nConds, nTRs_out);
                alld = nan(length(ROI_names), nConds,  nTRs_out);
                allconf = nan(length(ROI_names), nTrialsTotal, nTRs_out);
                allacc_rand = nan(length(ROI_names), nConds, nTRs_out, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nTRs_out, nPermIter);
            end

            %% voxel selection from each training set 
            % usng the training set only, find voxels that are selective
            % for response (which finger?)
            if ~isempty(nVox2Use) && nVox2Use<size(trnDat,2)
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                
                pvals = zeros(size(trnDat,2), 1);
                dat = trnDat;
                lab = trnLabs;
                parfor vx = 1:size(trnDat,2)
                     % choose the voxels        
                   [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                   pvals(vx) = pvalue;
                end
                assert(~any(isnan(pvals)))
                % use the correct column of pTable to sort
                [~, ind] = sort(pvals, 'ascend'); 
                voxelindsuse = ind(1:nVox2Use);
            else            
                % put in a placeholder here because using all voxels
                voxelindsuse = 1:size(trnDat,2);
            end

            for tr = 1:nTRs_out


                %% run the classifier

                %take out just data from this TR of interest.
                tstDat = squeeze(mainDat(:,tr,:));
               
                [predLabs,normEucDist] = normEucDistClass(trnDat(:,voxelindsuse),...
                    squeeze(tstDat(:,voxelindsuse)),trnLabs);

                acc = mean(predLabs==tstLabs);
                dprime = get_dprime(predLabs, tstLabs,tstLabs);

                allacc(vv,cc,tr) = acc;
                alld(vv,cc,tr) = dprime;
                
                % confidence is the distance to incorrect - distance to
                % correct. want a positive number (far from incorrect)
                conf = normEucDist(:,2) - normEucDist(:,1);
                conf(tstLabs==2) = -conf(tstLabs==2);
                % check these confidence labels to make sure they track -
                % always positive when classifier is correct, negative when
                % classifier makes a mistake.            
                assert(all(conf(predLabs==tstLabs)>0) && all(conf(predLabs~=tstLabs)<0))
                allconf(vv,trials2use,tr) = conf;

                if ismember(vv,areas2test)
                    fprintf('%s %s cc=%d tr=%d, performance on real data is %.2f, starting random shuffles over %d iters...\n',...
                    substr,ROI_names{vv},cc,tr,allacc(vv,cc,tr),nPermIter)

                    % now doing the permutation test, shuffle labels 1000 times.
                    randaccs= nan(nPermIter, 1);              
                    randd = nan(nPermIter, 1);
                    
                    % doing the shuffling before parfor loop
                    randlabs_all = zeros(size(trnLabs,1),nPermIter);
                    for ii=1:nPermIter
                        randlabs_all(:,ii) = trnLabs(randperm(numel(trnLabs)));                       
                    end
                    parfor ii=1:nPermIter

                        randlabs=randlabs_all(:,ii);
                            
                        [predLabs,normEucDist] = normEucDistClass(trnDat(:,voxelindsuse),...
                            squeeze(tstDat(:,voxelindsuse)), randlabs);

                        % get performance in each condition, for the random decoder
                        randaccs(ii) = mean(predLabs==tstLabs);                  
                        randd(ii) = get_dprime(predLabs,tstLabs,unique(tstLabs));

                    end
                   
                else
                    fprintf('%s %s cc=%d tr=%d, performance on real data is %.2f, skipping permutation test...\n',...
                    substr,ROI_names{vv},cc,tr,allacc(vv,cc,tr))
                    randaccs= nan(nPermIter, 1);              
                    randd = nan(nPermIter, 1);
                end
                
                % put everything into a big array for saving
                allacc_rand(vv,cc,tr,:) = randaccs;               
                alld_rand(vv,cc,tr,:) = randd;
                
            end
        end

    end

    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allacc_rand','alld_rand','allconf','condlabs','rt','correct');

end