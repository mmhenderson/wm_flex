%% Response decoding analysis
% Train and test linear decoder, using data from one task condition at a
% time. Labels are the expected (correct) response on each trial.
% Train/test within and across timepts, giving full [nTRs x nTRs] matrix
% cross-validate across sessions, response/luminance mapping was swapped
% between sessions.
% Saves the results in a mat file, can plot it using a separate script
% (plotClassResults_Response_AcrossTime.m)
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

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

nVox2Use = 10000;
nPermIter=1000;

class_str = 'normEucDist';

nTrialsTotal=400;

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponse_AcrossTime_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

    for vv = 1:length(mainSig)
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        %% loop over conditions
        for cc =1:nConds

            respLabs = mainSig(vv).CorrectResp;
            condLabs = mainSig(vv).condLabs;
            runLabs = mainSig(vv).runLabs;
            % take out just relevant condition here
            trials2use = ~isnan(respLabs) & respLabs~=0 & condLabs==cc;

            respLabs = respLabs(trials2use);
            runLabs = runLabs(trials2use);

            sessLabs = ones(size(runLabs));
            sessLabs(runLabs>10) = 2;
            cvLabs = sessLabs;
            
            mainDat = mainSig(vv).dat_by_TR;
            mainDat = mainDat(trials2use,:,:);

            % mainDat is [ntrials x nTRs x nVox]
            nTRs_out = size(mainDat,2);

            if vv==1 && cc==1
                % preallocate array here
                allacc = nan(length(ROI_names), nConds, nTRs_out, nTRs_out);
                alld = nan(length(ROI_names), nConds,  nTRs_out, nTRs_out);
                allconf = nan(length(ROI_names), nTrialsTotal, nTRs_out, nTRs_out);
                allacc_rand = nan(length(ROI_names), nConds, nTRs_out, nTRs_out, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nTRs_out, nTRs_out, nPermIter);
            end

            nCV = numel(unique(cvLabs));

            for tr1 = 1:nTRs_out
               
                % take out just data from this TR of interest.
                dat2use_trn = squeeze(mainDat(:,tr1,:));

                %% voxel selection from each training set 
                % for this voxel selection I'm using trials from all conditions, but
                % leaving out one session at a time. This gives a list of voxels to use
                % for each fold of cross validation. The same voxels are used
                % regardless of which condition we are using for classification. Think
                % this makes the condition comparisons more fair. Also saves time
                % because we only need to run this once.
                if ~isempty(nVox2Use) && nVox2Use<size(dat2use_trn,2)
                    % get ready for parallel pool operations
                    numcores = 8;
                    if isempty(gcp('nocreate'))
                        parpool(numcores);
                    end
                    fprintf('running voxel selection f-test for %s %s: %s condition, tr=%d\n',substr, ROI_names{vv},condLabStrs{cc}, tr1)
                    voxStatTable = zeros(size(dat2use_trn,2),nCV);
                    for rr = 1:nCV
                        inds = cvLabs~=rr;
                        pvals = zeros(size(dat2use_trn,2), 1);
                        dat = dat2use_trn(inds,:);
                        lab = respLabs(inds,:);
                        parfor vx = 1:size(dat2use_trn,2)
                             % choose the voxels        
                           [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                           pvals(vx) = pvalue;
                        end 
                        voxStatTable(:,rr) = pvals;
                    end
                    nVox2Use_now = nVox2Use;
                else            
                    % put in a placeholder here because using all voxels
                    voxStatTable = zeros(size(dat2use_trn,2),nCV);
                    nVox2Use_now = [];
                end

                for tr2 = 1:nTRs_out
                    
                    dat2use_tst = squeeze(mainDat(:,tr2,:));
                    
                    %% run the classifier 

                    trnDat = dat2use_trn;
                    trnLabs = respLabs;
                    trnCV = cvLabs;

                    tstDat = dat2use_tst;
                    tstLabs = respLabs;
                    tstCV = cvLabs;
                    
                    % using custom code to do cross-validation 
                    [~,~,predLabs,normEucDist] = my_classifier_cross_wconf(trnDat,trnLabs,...
                        trnCV,tstDat, tstLabs,...
                        tstCV,class_str,100,nVox2Use_now,voxStatTable,1);

                    acc = mean(predLabs==tstLabs);
                    dprime = get_dprime(predLabs, tstLabs,tstLabs);

                    allacc(vv,cc,tr1,tr2) = acc;
                    alld(vv,cc,tr1,tr2) = dprime;

                    % confidence is the distance to incorrect - distance to
                    % correct. want a positive number (far from incorrect)
                    conf = normEucDist(:,2) - normEucDist(:,1);
                    conf(tstLabs==2) = -conf(tstLabs==2);
                    % check these confidence labels to make sure they track -
                    % always positive when classifier is correct, negative when
                    % classifier makes a mistake.            
                    assert(all(conf(predLabs==tstLabs)>0) && all(conf(predLabs~=tstLabs)<0))
                    allconf(vv,trials2use,tr1,tr2) = conf;

                end
            end
        end

    end

    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allacc_rand','alld_rand','allconf','condlabs','rt','correct');

end