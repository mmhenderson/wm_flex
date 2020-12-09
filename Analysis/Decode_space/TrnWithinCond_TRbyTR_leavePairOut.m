% MMH 10/29/20

% training/testing on SPATIAL WORKING MEMORY localizer, testing on task data.

clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
addpath(fullfile(exp_path,'Analysis','stats_code'));

nVox2Use = 10000;    % this is the max number of vox to use, so if it's very big we're using all the voxels

nPermIter=1000;

class_str = 'normEucDist';
% class_str = 'svmtrain_lin';

% dbstop if error
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 456566;
rng(rndseed,'twister');

nBins=8;
bin_centers=[0:45:359];
bin_size=diff(bin_centers(1:2));
% going to do 4 different classifications - bin 1 versus bin 5, 2
% versus 6, etc.
groups = [(1:4)',(5:8)'];
nGroups = size(groups,1);


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

ips_inds = [6:9];
ips01_inds=[6,7];
ips23_inds=[8,9];
for ss=1:length(sublist)
    
    allchanresp = [];

    substr = sprintf('S%02d',sublist(ss));
%     fn2load = fullfile(exp_path,'Samples',sprintf('SWMLocSignalByTrial_%s.mat',substr));
%     load(fn2load);
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnWithinCond_TRbyTR_leavePairOut_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
    %% loop over ROIs and run the model for each.
    ROI_names{end+1} = 'IPS0-3';
    ROI_names{end+1} = 'IPS0-1';
    ROI_names{end+1} = 'IPS2-3';
    for vv = 1:length(ROI_names) 
        
        %% create the merged IPS regions
        if vv==length(ROI_names)
            mainSig(vv) = mainSig(ips23_inds(1));
            for ip = ips23_inds(2:end)
                mainSig(vv).dat_avg = cat(2, mainSig(vv).dat_avg, mainSig(ip).dat_avg);
                mainSig(vv).dat_by_TR = cat(3, mainSig(vv).dat_by_TR, mainSig(ip).dat_by_TR);
            end                
        end
        if vv==length(ROI_names)-1
            mainSig(vv) = mainSig(ips01_inds(1));
            for ip = ips01_inds(2:end)
                mainSig(vv).dat_avg = cat(2, mainSig(vv).dat_avg, mainSig(ip).dat_avg);
                mainSig(vv).dat_by_TR = cat(3, mainSig(vv).dat_by_TR, mainSig(ip).dat_by_TR);
            end                
        end
        if vv==length(ROI_names)-2
            mainSig(vv) = mainSig(ips_inds(1));
            for ip = ips_inds(2:end)
                mainSig(vv).dat_avg = cat(2, mainSig(vv).dat_avg, mainSig(ip).dat_avg);
                mainSig(vv).dat_by_TR = cat(3, mainSig(vv).dat_by_TR, mainSig(ip).dat_by_TR);
            end                
        end
         %% check if there's data here

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<2
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s\n',ROI_names{vv});
        
        for cc=1:nConds
            
            % taking out one condition at a time
            condLabs = mainSig(vv).condLabs;
            inds2use = condLabs==cc;

            posLabs = mainSig(vv).targPos(inds2use,:);
            allDat = mainSig(vv).dat_by_TR(inds2use,:,:); % [nTrials x nTRs x nVox]
            % trying to get rid of baseline shifts here. do subtraction within
            % a TR only. 
            allDat = allDat - repmat(mean(allDat,3), 1, 1, size(allDat,3));
    %         
            runLabs = mainSig(vv).runLabs;
            runLabs = runLabs(inds2use);

%             sessLabs = ones(size(runLabs));
%             sessLabs(runLabs>10) = 2;
%             cvLabs = sessLabs;
%             nCV = numel(unique(cvLabs));
            
            if vv==1 && cc==1
                % preallocate array here
                nTRs=size(allDat,2);
                allacc = nan(length(ROI_names), nConds, nGroups, nTRs);
                alld = nan(length(ROI_names), nConds, nGroups, nTRs);
                allconf = nan(length(ROI_names), numel(condLabs), nTRs);
                allacc_rand = nan(length(ROI_names), nConds, nGroups, nTRs, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nGroups, nTRs, nPermIter);
            end

            conf_this_cond = nan(numel(posLabs),nTRs);
            
            % bin these for classifier - want 8 bins that are roughly centered at
            % 0, 45, 90, 135,...
            binLabs = zeros(size(posLabs));        
            for bb=1:nBins
                inds_this_bin = abs(posLabs-(bin_centers(bb)-0.0001))<bin_size/2 | abs((posLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
                binLabs(inds_this_bin) = bb;
            end
            neach_tst = sum(repmat(binLabs,1,nBins)==repmat((1:nBins),size(binLabs,1),1));
            assert(~any(binLabs==0))

            % create cross-validation labels
            % on each fold, want to leave out one pair of trials (one from
            % each label bin) so that training set stays balanced.
            % basically each bin contributes once to each cross-validation
            % fold.           
            cvLabs=zeros(size(binLabs));
            for bb=1:nBins
                inds = binLabs==bb;
                cvLabs(inds) = 1:sum(inds);
            end
            assert(~any(cvLabs==0))
            nCV = numel(unique(cvLabs));
            
            for tr=1:nTRs

                % take out just data from this TR of interest.
                dat2use = squeeze(allDat(:,tr,:));

                %% voxel selection from the training set 

                if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                    fprintf('running voxel selection f-test for %s %s: %s condition, tr=%d\n',substr, ROI_names{vv},condLabStrs{cc}, tr)
                    voxStatTable = zeros(size(dat2use,2),nCV);
                    for rr = 1:nCV
                        inds = cvLabs~=rr;
                        pvals = zeros(size(dat2use,2), 1);
                        dat = dat2use(inds,:);
                        lab = binLabs(inds,:);
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

                %% run each classification scheme

                for gg=1:nGroups

                    % taking out only trn and test trials with the relevant two
                    % labels that we're trying to discriminate
                    inds2use = ismember(binLabs, groups(gg,:));
                
                    group_labs = binLabs(inds2use);
                    dat_this_group = squeeze(dat2use(inds2use, :));
                    cv_this_group = cvLabs(inds2use);
                    
                    runs_using = unique(cv_this_group);
                    % train/test the decoder
                    [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(dat_this_group,group_labs,...
                        cv_this_group,dat_this_group, group_labs,...
                        cv_this_group,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),1);
%                     [thesepredlabs,normEucDist] = normEucDistClass(allDatUse,tstDatUse,reallabs);

                    % compute accuracy on the subset of test trials in each
                    % condition separately                  
                    allacc(vv,cc,gg,tr) = mean(thesepredlabs==group_labs);
                    alld(vv,cc,gg,tr) = get_dprime(thesepredlabs,group_labs,unique(group_labs));
                  
                    % confidence is the distance to incorrect - distance to
                    % correct. want a positive number (far from incorrect)
                    conf = normEucDist(:,2) - normEucDist(:,1);
                    conf(group_labs==groups(gg,2)) = -conf(group_labs==groups(gg,2));
                    % check these confidence labels to make sure they track -
                    % always positive when classifier is correct, negative when
                    % classifier makes a mistake.
                    assert(all(conf(thesepredlabs==group_labs)>0) && all(conf(thesepredlabs~=group_labs)<0))
                    conf_this_cond(inds2use,tr) = conf;

                    fprintf('%s %s %s: gg=%d, tr=%d, performance on real data is %.2f, starting random shuffles over %d iters...\n',...
                        substr,ROI_names{vv},condLabStrs{cc},gg,tr,allacc(vv,cc,gg,tr),nPermIter)

                    % now doing the permutation test, shuffle labels 1000 times.
                    randaccs = nan(nPermIter, 1);
                    randd = nan(nPermIter, 1);
                  
                    parfor ii=1:nPermIter
                        % randomize all labels (note this is across all runs,
                        % so we're shuffling training and testing sets at once.
                        randlabs=nan(size(group_labs));
                        for cv=1:numel(runs_using)
                            % shuffle the data from one cross-validation fold at a time, so we
                            % don't un-balance the training sets. 
                            inds=cv_this_group==runs_using(cv);
                            dat2shuff=group_labs(inds);
                            randlabs(inds) = dat2shuff(randperm(numel(dat2shuff)));
                        end                  
                        % run classifier with the random labels
                        % train/test the decoder
                        [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(dat_this_group,randlabs,...
                            cv_this_group,dat_this_group, randlabs,...
                            cv_this_group,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),1);
                        % get performance in each condition, for the random decoder
                        randaccs(ii) = mean(thesepredlabs==group_labs);                        
                        randd(ii) = get_dprime(thesepredlabs,group_labs,unique(group_labs));
                       
                    end

                    % put everything into a big array for saving
                    allacc_rand(vv,cc,gg,tr,:) = randaccs;                   
                    alld_rand(vv,cc,gg,tr,:) = randd;
                  
                end
            end
            
            % put confidence labels back together over all conditions
            allconf(vv, condLabs==cc,:) = conf_this_cond;
        end
    end
    
    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','rt','correct','condlabs','allacc_rand','alld_rand');

end