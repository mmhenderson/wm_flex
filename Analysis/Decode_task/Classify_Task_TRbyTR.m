%% Task decoding analysis
% Train and test linear decoder, using data from all trials of main task.
% Labels are the conditions - random and predictable.
% Train and test within TR for time-resolved decoding.
% Saves the results in a mat file, can plot it using a separate script
% (plotClassResults_Task_TRbyTR.m)
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
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 686787;
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
    fn2save = fullfile(save_dir,sprintf('ClassifyTask_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

    areas2test = [1:length(mainSig)];
   
    for vv = 1:length(mainSig)

        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        mainDat = mainSig(vv).dat_by_TR;
        % mainDat is [ntrials x nTRs x nVox]
        nTRs_out = size(mainDat,2);
        % subtract mean over voxels (third dim)
        mainDat = mainDat - repmat(mean(mainDat,3), 1, 1, size(mainDat, 3));
 
        if vv==1
            % preallocate array here
            allacc = nan(length(ROI_names), nTRs_out);
            alld = nan(length(ROI_names), nTRs_out);
            allacc_rand = nan(length(ROI_names), nTRs_out, nPermIter);
            alld_rand = nan(length(ROI_names), nTRs_out, nPermIter);
        end
        
        condLabs = mainSig(vv).condLabs;
        runLabs = mainSig(vv).runLabs;
        
        assert(numel(unique(runLabs))==20);
        sessLabs = ones(size(runLabs));
        sessLabs(runLabs>10) = 2;
        cvLabs = sessLabs;
        nCV = numel(unique(cvLabs));
        
        for tr = 1:nTRs_out
            
            % take out just data from this TR of interest.
            dat2use = squeeze(mainDat(:,tr,:));
            
            %% voxel selection from each training set 
            % for this voxel selection I'm using trials from all conditions, but
            % leaving out one session at a time. This gives a list of voxels to use
            % for each fold of cross validation. The same voxels are used
            % regardless of which condition we are using for classification. Think
            % this makes the condition comparisons more fair. Also saves time
            % because we only need to run this once.
            if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                fprintf('running voxel selection f-test for %s %s, tr=%d\n',substr, ROI_names{vv},tr)
                voxStatTable = zeros(size(dat2use,2),nCV);
                for rr = 1:nCV
                    inds = cvLabs~=rr;
                    pvals = zeros(size(dat2use,2), 1);
                    dat = dat2use(inds,:);
                    lab = condLabs(inds,:);
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

            %% run classifier
            trnDat = dat2use;
            trnLabs = condLabs;
            trnCV = cvLabs;

            tstDat = dat2use;
            tstLabs = condLabs;
            tstCV = cvLabs;
            
            % train/test the decoder - use a custom function to do
            % cross-validation here. 
            [~,~,predLabs] = my_classifier_cross(trnDat,trnLabs,...
                trnCV,tstDat, tstLabs,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,0);

            acc = mean(predLabs==tstLabs);
            dprime = get_dprime(predLabs, tstLabs,tstLabs);

            allacc(vv,tr) = acc;
            alld(vv,tr) = dprime;
            
    
            if ismember(vv,areas2test)
                fprintf('%s %s tr=%d, performance on real data is %.2f, starting random shuffles over %d iters...\n',...
                substr,ROI_names{vv},tr,allacc(vv,tr),nPermIter)

                % now doing the permutation test, shuffle labels 1000 times.
                randaccs= nan(nPermIter, 1);              
                randd = nan(nPermIter, 1);

                % doing the shuffling before parfor loop 
                randlabs_all = zeros(size(trnLabs,1),nPermIter);
                for ii=1:nPermIter
                     for se=1:2
                        % shuffle the data from one session at a time, so we
                        % don't un-balance the training sets. 
                        inds=trnCV==se;
                        dat2shuff=trnLabs(inds);
                        randlabs_all(inds,ii) = dat2shuff(randperm(numel(dat2shuff)));
                    end                 
                end  
                parfor ii=1:nPermIter
                   randlabs=randlabs_all(:,ii);          
                    % run classifier with the random labels
                    [~,~,predLabs] = my_classifier_cross(trnDat,randlabs,...
                    trnCV,tstDat, randlabs,...
                    tstCV,class_str,100,nVox2Use_now,voxStatTable,0);

                    % get performance in each condition, for the random decoder
                    randaccs(ii) = mean(predLabs==randlabs);                  
                    randd(ii) = get_dprime(predLabs,randlabs,unique(randlabs));

                end
       
            else
                fprintf('%s %s tr=%d, performance on real data is %.2f, skipping permutation test...\n',...
                substr,ROI_names{vv},tr,allacc(vv,tr))
                randaccs= nan(nPermIter, 1);              
                randd = nan(nPermIter, 1);
            end

            % put everything into a big array for saving
            allacc_rand(vv,tr,:) = randaccs;               
            alld_rand(vv,tr,:) = randd;

            
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allacc_rand','alld_rand');

end