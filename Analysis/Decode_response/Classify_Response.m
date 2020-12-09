% MMH 3/12/20
% classifying RESPONSE - WHICH FINGER DID THEY EVENTUALLY PRESS??
    % note this is counter-balanced w/r/t COLOR of the side.

clear
close all;

sublist = [7];
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;
filesepinds = find(curr_dir==filesep);
nDirsUp = 2;
exp_path = curr_dir(1:filesepinds(end-nDirsUp+1));
addpath(fullfile(exp_path,'Analysis','stats_code'));

nVox2Use = 10000;
nPermIter=1000;
condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

class_str = 'normEucDist';
% class_str = 'svmtrain_lin';
nTrialsTotal=2*10*20;

dbstop if error
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end

rndseed = 346456;
rng(rndseed,'twister');
%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

   
    v2do = 1:length(ROI_names);
    
    allconf = nan(length(ROI_names), nTrialsTotal);
    
    for vv = v2do
%     for vv = 1:length(ROI_names)
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        for cc = 1:nConds
            
%             respLabs = mainSig(vv).RespActual;
            respLabs = mainSig(vv).CorrectResp;
            condLabs = mainSig(vv).condLabs;
            runLabs = mainSig(vv).runLabs;
            % getting rid of any trials with no response here
            % also taking out just the relevant condition!!
            trials2use = ~isnan(respLabs) & respLabs~=0 & condLabs==cc;
            
            respLabs = respLabs(trials2use);
            runLabs = runLabs(trials2use);
            sessLabs = ones(size(runLabs));
            sessLabs(runLabs>10) = 2;
            cvLabs = sessLabs;
            
            mainDat = mainSig(vv).dat_avg;
            mainDat = mainDat(trials2use,:);

            assert(sum(respLabs==1 & sessLabs==1)==50);
            assert(sum(respLabs==1 & sessLabs==2)==50);
            % subtract mean over voxels 
%             mainDat = mainDat - repmat(mean(mainDat,2), 1, size(mainDat, 2));

            if vv==v2do(1) && cc==1
                % preallocate array here
                allacc = nan(length(ROI_names), nConds);
                alld = nan(length(ROI_names), nConds);
                allacc_rand = nan(length(ROI_names), nConds, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nPermIter);
            end

            nCV = numel(unique(cvLabs));

            dat2use = mainDat;

            %% voxel selection from each training set 
            % for this voxel selection I'm using trials from all conditions, but
            % leaving out one session at a time. This gives a list of voxels to use
            % for each fold of cross validation. The same voxels are used
            % regardless of which condition we are using for classification. Think
            % this makes the condition comparisons more fair. Also saves time
            % because we only need to run this once.
            if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                voxStatTable = zeros(size(dat2use,2),nCV);
                for rr = 1:nCV
                    inds = cvLabs~=rr;
                    pvals = zeros(size(dat2use,2), 1);
                    dat = dat2use(inds,:);
                    lab = respLabs(inds,:);
                    parfor vx = 1:size(dat2use,2)
                         % choose the voxels        
                       [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                       pvals(vx) = pvalue;
                    end 
                    voxStatTable(:,rr) = pvals;
                end
                nVox2Use_now = nVox2Use;
            else            
                % put in a placeholder here because using all voxels
                voxStatTable = zeros(size(dat2use,2),nCV);
                nVox2Use_now = [];
            end

            %% define train and test set 

            % same data here because we're not cross-generalizing or anything
            trnDat = dat2use;
            trnLabs = respLabs;
            trnCV = cvLabs;

            tstDat = dat2use;
            tstLabs = respLabs;
            tstCV = cvLabs;

            %% run the classifier w/ balancing if needed
            
            [~,~,predLabs,normEucDist] = my_classifier_cross_wconf(trnDat,trnLabs,...
                trnCV,tstDat, tstLabs,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,1);

            acc = mean(predLabs==tstLabs);
            dprime = get_dprime(predLabs, tstLabs,tstLabs);
            
            % confidence is the distance to incorrect - distance to
            % correct. want a positive number (far from incorrect)
            conf = normEucDist(:,2) - normEucDist(:,1);
            conf(tstLabs==2) = -conf(tstLabs==2);
            % check these confidence labels to make sure they track -
            % always positive when classifier is correct, negative when
            % classifier makes a mistake.            
            assert(all(conf(predLabs==tstLabs)>0) && all(conf(predLabs~=tstLabs)<0))
            allconf(vv,trials2use) = conf;
            
            allacc(vv,cc) = acc;
            alld(vv,cc) = dprime;
            
            %% now doing the permutation test, shuffle labels 1000 times.
            randaccs= nan(nPermIter, 1);              
            randd = nan(nPermIter, 1);

            parfor ii=1:nPermIter
                % randomize all labels (note this is across all runs,
                % so we're shuffling training and testing sets at once.
                randlabs_all=nan(size(trnLabs));
                for se=1:2
                    % shuffle the data from one session at a time, so we
                    % don't un-balance the training sets. 
                    inds=trnCV==se;
                    dat2shuff=trnLabs(inds);
                    randlabs_all(inds) = dat2shuff(randperm(numel(dat2shuff)));
                end
                % run classifier with the random labels
                [~,~,predLabs] = my_classifier_cross(trnDat,randlabs_all,...
                trnCV,tstDat, randlabs_all,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,1);

                % get performance in each condition, for the random decoder
                randaccs(ii) = mean(predLabs==randlabs_all);                  
                randd(ii) = get_dprime(predLabs,randlabs_all,unique(randlabs_all));

            end

            % put everything into a big array for saving
            allacc_rand(vv,cc,:) = randaccs;               
            alld_rand(vv,cc,:) = randd;
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allacc_rand','alld_rand','allconf');

end