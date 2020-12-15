%% Spatial decoding analysis
% Train and test linear decoder,using data from one task at a time (leave
% out a small number of trials at a time).Train/test at each TR for time 
% resolved decoding accuracy.
% grouping spatial positions into 8 bins 45 deg wide, then doing binary
% decoding between bins that are 180 deg apart (e.g. 0 vs 180). 
% Saves the results in a mat file, can plot it using a separate script
% (plotTrnWithinCond_TRbyTR.m)

%%
clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
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
rndseed = 456566;
rng(rndseed,'twister');

% define the spatial position bins
nBins=8;
bin_centers=[0:45:359];
bin_size=diff(bin_centers(1:2));
% going to do 4 different classifications - bin 1 versus bin 5, 2
% versus 6, etc.
groups = [(1:4)',(5:8)'];
nGroups = size(groups,1);


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

% also going to do this analyis for "merged" IPS ROIs that combine multiple
% subregions of IPS at a time. 
% these get done after all the individual ROIs.
ips_inds = [6:9];
ips01_inds=[6,7];
ips23_inds=[8,9];
for ss=1:length(sublist)
    
    allchanresp = [];

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnWithinCond_TRbyTR_leavePairOut_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
    ROI_names{end+1} = 'IPS0-3';
    ROI_names{end+1} = 'IPS0-1';
    ROI_names{end+1} = 'IPS2-3';
    
    %% loop over ROIs and run the model for each.
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
            % subtract mean over voxels. do subtraction within
            % a TR only. 
            allDat = allDat - repmat(mean(allDat,3), 1, 1, size(allDat,3));

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
            % note that the number of trials in each of these 8 bins is not
            % exactly equal - however, the bins that are 180 deg apart
            % always have the same number of trials. So each decoder always
            % has perfectly balanced training set.
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
                    % train/test the decoder - use a custom function to do
                    % cross-validation here. 
                    [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(dat_this_group,group_labs,...
                        cv_this_group,dat_this_group, group_labs,...
                        cv_this_group,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),0);

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
                    % doing the shuffling before parfor loop
                    randlabs_all = zeros(size(group_labs,1),nPermIter);
                    for ii=1:nPermIter
                        for cv=1:numel(runs_using)
                            % shuffle the data from one cross-validation fold at a time, so we
                            % don't un-balance the training sets.
                            inds=cv_this_group==runs_using(cv);
                            dat2shuff=group_labs(inds);
                            randlabs_all(inds,ii) = dat2shuff(randperm(numel(dat2shuff)));
                        end
                    end
                    parfor ii=1:nPermIter
                        randlabs=randlabs_all(:,ii);
                        % run classifier with the random labels
                        [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(dat_this_group,randlabs,...
                            cv_this_group,dat_this_group, randlabs,...
                            cv_this_group,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),0);
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