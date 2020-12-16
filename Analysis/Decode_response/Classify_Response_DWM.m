%% Response decoding analysis
% Train and test linear decoder, using data from digit working memory (DWM)
% task (button pressing with delayed response).
% Labels are the expected (correct) response on each trial.
% cross-validate across sessions, response/luminance mapping was swapped
% between sessions.
% Saves the results in a mat file, can plot it using a separate script
% (plotClassResults_ResponseDWMLoc.m)
% Note this task isn't included in our paper
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
rndseed = 121343;
rng(rndseed,'twister');

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);
nTrialsTotal=2*10*20;

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('DWMLocSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponse_TrnTestDWMLoc_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

   
    v2do = 1:length(ROI_names);
    
    
    nTrialsTotal = size(locSig(1).dat_avg_targ,1);
    allconf = nan(length(ROI_names), nTrialsTotal);
    for vv = v2do

        %% pull out the data for main task

        if length(locSig)<vv || isempty(locSig(vv).dat_avg_targ) || size(locSig(vv).dat_avg_targ,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        for cc = 1:nConds
            
            % using all trials because accuracy is generally very high. Can
            % also remove the incorrect trials, but that would un-balance
            % the training set.
            trials2use = ones(size(locSig(1).ActDigit,1),1)==1;
            respLabs = locSig(vv).ExpDigit(trials2use);
            runLabs = locSig(vv).runLabs(trials2use);
            digDat = locSig(vv).dat_avg_targ(trials2use,:);

            if vv==v2do(1)
                % preallocate array here
                allacc = nan(length(ROI_names), nConds);
                alld = nan(length(ROI_names), nConds);
                allacc_rand = nan(length(ROI_names), nConds, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nPermIter);
            end

            nRuns = numel(unique(runLabs));

            dat2use = digDat;

            %% voxel selection from each training set 
            % for this voxel selection I'm using trials from all conditions, but
            % leaving out one session at a time. This gives a list of voxels to use
            % for each fold of cross validation. The same voxels are used
            % regardless of which condition we are using for classification. Think
            % this makes the condition comparisons more fair. Also saves time
            % because we only need to run this once.
            if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                voxStatTable = zeros(size(dat2use,2),nRuns);
                for rr = 1:nRuns
                    inds = runLabs~=rr;
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
                voxStatTable = zeros(size(dat2use,2),nRuns);
                nVox2Use_now = [];
            end

            %% run the classifier 
            
            trnDat = dat2use;
            trnLabs = respLabs;
            trnCV = runLabs;

            tstDat = dat2use;
            tstLabs = respLabs;
            tstCV = runLabs;
            
            % using custom code to do cross-validation - same data goes in
            % as train and test, but cross-validation labels determines
            % which part used to train and test.
            [~,~,predLabs,normEucDist] = my_classifier_cross_wconf(trnDat,trnLabs,...
                trnCV,tstDat, tstLabs,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,0);

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

            % doing the shuffling before parfor loop 
            randlabs_all = zeros(size(trnLabs,1),nPermIter);
            for ii=1:nPermIter
                 for se=1:nRuns
                    % shuffle the data from one session at a time, so we
                    % don't un-balance the training sets. 
                    inds=trnCV==se;
                    dat2shuff=trnLabs(inds);
                    randlabs_all(inds,ii) = dat2shuff(randperm(numel(dat2shuff)));
                 end                
            end  
            parfor ii=1:nPermIter 
                randlabs=randlabs_all(:,ii)                
                % run classifier with the random labels
                [~,~,predLabs] = my_classifier_cross(trnDat,randlabs,...
                trnCV,tstDat, randlabs,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,0);
                % get performance in each condition, for the random decoder
                randaccs(ii) = mean(predLabs==randlabs);                  
                randd(ii) = get_dprime(predLabs,randlabs,unique(randlabs));
            end
           
            % put everything into a big array for saving
            allacc_rand(vv,cc,:) = randaccs;               
            alld_rand(vv,cc,:) = randd;
            
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','allacc_rand','alld_rand');

end