%% Spatial decoding analysis
% Train and test linear decoder, using data from one task at a time (leave
% out a small number of trials at a time). Train and test across timepts
% (each TR serves as a training set for each other TR to give a full [nTR x
% nTR] matrix)
% grouping spatial positions into 8 bins 45 deg wide, then doing binary
% decoding between bins that are 180 deg apart (e.g. 0 vs 180). 
% Saves the results in a mat file, can plot it using a separate script
% (plotTrnWithinCond_AcrossTime.m)
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

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnWithinCond_AcrossTime_leavePairOut_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
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
            % remove mean over voxels, within each TR
            allDat = allDat - repmat(mean(allDat,3), 1, 1, size(allDat,3));

            if vv==1 && cc==1
                % preallocate array here
                nTRs=size(allDat,2);
                allacc = nan(length(ROI_names), nConds, nGroups, nTRs, nTRs);
                alld = nan(length(ROI_names), nConds, nGroups, nTRs, nTRs);
            end

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
            
            %% loop over trs to train model
            for tr_trn=1:nTRs
                
                % take out just data from this TR of interest.
                dat2use_trn = squeeze(allDat(:,tr_trn,:));
                
                %% voxel selection from the training set 
                if ~isempty(nVox2Use) && nVox2Use<size(dat2use_trn,2)
                    % get ready for parallel pool operations
                    numcores = 8;
                    if isempty(gcp('nocreate'))
                        parpool(numcores);
                    end
                    fprintf('running voxel selection f-test for %s %s: %s condition, train tr=%d\n',substr, ROI_names{vv},condLabStrs{cc}, tr_trn)
                    voxStatTable = zeros(size(dat2use_trn,2),nCV);
                    for rr = 1:nCV
                        inds = cvLabs~=rr;
                        pvals = zeros(size(dat2use_trn,2), 1);
                        dat = dat2use_trn(inds,:);
                        lab = binLabs(inds,:);
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
                
                %% loop over TRs to test model
                for tr_tst = 1:nTRs 

                    % take out just data from this TR of interest.
                    dat2use_tst = squeeze(allDat(:,tr_tst,:));

                    %% run each classification scheme

                    for gg=1:nGroups

                        fprintf('%s %s %s, train tr=%d test tr=%d, group %d\n',substr,ROI_names{vv}, condLabStrs{cc},tr_trn, tr_tst, gg)
                        % taking out only trn and test trials with the relevant two
                        % labels that we're trying to discriminate
                        inds2use = ismember(binLabs, groups(gg,:));

                        group_labs = binLabs(inds2use);
                        trndat_this_group = squeeze(dat2use_trn(inds2use, :));
                        tstdat_this_group = squeeze(dat2use_tst(inds2use, :));
                        cv_this_group = cvLabs(inds2use);

                        runs_using = unique(cv_this_group);
                        % train/test the decoder
                        [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(trndat_this_group,group_labs,...
                            cv_this_group,tstdat_this_group, group_labs,...
                            cv_this_group,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),1);

                        % compute accuracy                 
                        allacc(vv,cc,gg,tr_trn,tr_tst) = mean(thesepredlabs==group_labs);
                        alld(vv,cc,gg,tr_trn,tr_tst) = get_dprime(thesepredlabs,group_labs,unique(group_labs));


                    end
                end
            end
           
        end
    end
    
    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','rt','correct','condlabs');

end