%% Spatial decoding analysis
% Train and test linear decoder, across conditions - so one task served as
% training and one task as testing set. 
% grouping spatial positions into 8 bins 45 deg wide, then doing binary
% decoding between bins that are 180 deg apart (e.g. 0 vs 180). 
% Saves the results in a mat file, can plot it using a separate script
% (plotTrnTestWithinAcrossCond.m)
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
    fn2save = fullfile(save_dir,sprintf('TrnAcrossCond_leavePairOut_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
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
            
            % test and train are OPPOSITE conditions
            condLabs = mainSig(vv).condLabs;
            inds2use_tst = condLabs==cc;    
            inds2use_trn = condLabs~=cc;

            posLabs_tst = mainSig(vv).targPos(inds2use_tst,:);
            posLabs_trn = mainSig(vv).targPos(inds2use_trn,:);
            allDat_tst = mainSig(vv).dat_avg(inds2use_tst,:); % [nTrials x nVox]
            allDat_trn = mainSig(vv).dat_avg(inds2use_trn,:); % [nTrials x nVox]
            % subtract the mean over voxels
            allDat_tst = allDat_tst - repmat(mean(allDat_tst,2), 1, size(allDat_tst,2));
            allDat_trn = allDat_trn - repmat(mean(allDat_trn,2), 1, size(allDat_trn,2));

            if vv==1 && cc==1
                % preallocate array here
                allacc = nan(length(ROI_names), nConds, nGroups);
                alld = nan(length(ROI_names), nConds, nGroups);
                allconf = nan(length(ROI_names), numel(condLabs));
                allacc_rand = nan(length(ROI_names), nConds, nGroups, nPermIter);
                alld_rand = nan(length(ROI_names), nConds, nGroups, nPermIter);
            end

            conf_this_cond = nan(numel(posLabs_tst),1);
            
            % bin these for classifier - want 8 bins that are roughly centered at
            % 0, 45, 90, 135,...
            binLabs_tst = zeros(size(posLabs_tst));        
            for bb=1:nBins
                inds_this_bin = abs(posLabs_tst-(bin_centers(bb)-0.0001))<bin_size/2 | abs((posLabs_tst-360)-(bin_centers(bb)-0.0001))<bin_size/2;
                binLabs_tst(inds_this_bin) = bb;
            end
            neach_tst = sum(repmat(binLabs_tst,1,nBins)==repmat((1:nBins),size(binLabs_tst,1),1));
            assert(~any(binLabs_tst==0))

            binLabs_trn = zeros(size(posLabs_trn));        
            for bb=1:nBins
                inds_this_bin = abs(posLabs_trn-(bin_centers(bb)-0.0001))<bin_size/2 | abs((posLabs_trn-360)-(bin_centers(bb)-0.0001))<bin_size/2;
                binLabs_trn(inds_this_bin) = bb;
            end
            neach_trn = sum(repmat(binLabs_trn,1,nBins)==repmat((1:nBins),size(binLabs_trn,1),1));
            assert(~any(binLabs_trn==0))

            % create cross-validation labels
            % on each fold, want to leave out one pair of trials (one from
            % each label bin) so that training set stays balanced.
            % basically each bin contributes once to each cross-validation
            % fold.           
            cvLabs_tst=zeros(size(binLabs_tst));
            for bb=1:nBins
                inds = binLabs_tst==bb;
                cvLabs_tst(inds) = 1:sum(inds);
            end
            assert(~any(cvLabs_tst==0))
            nCV = numel(unique(cvLabs_tst));
            
            cvLabs_trn=zeros(size(binLabs_trn));
            for bb=1:nBins
                inds = binLabs_trn==bb;
                cvLabs_trn(inds) = 1:sum(inds);
            end
            assert(~any(cvLabs_trn==0))
            assert(nCV==numel(unique(cvLabs_trn)))
            %% voxel selection from the training set 

            if ~isempty(nVox2Use) && nVox2Use<size(allDat_trn,2)
                % get ready for parallel pool operations
                numcores = 8;
                if isempty(gcp('nocreate'))
                    parpool(numcores);
                end
                fprintf('running voxel selection f-test for %s %s: %s condition\n',substr, ROI_names{vv},condLabStrs{cc})
                voxStatTable = zeros(size(allDat_trn,2),nCV);
                for rr = 1:nCV
                    inds = cvLabs_trn~=rr;
                    pvals = zeros(size(allDat_trn,2), 1);
                    dat = allDat_trn(inds,:);
                    lab = binLabs_trn(inds,:);
                    parfor vx = 1:size(allDat_trn,2)
                         % choose the voxels        
                       [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                       pvals(vx) = pvalue;
                    end 
                    voxStatTable(:,rr) = pvals;
                end
                nVox2Use_now = nVox2Use;
            else            
                % put in a placeholder here because using all voxels
                voxStatTable = zeros(size(allDat_trn,2),nCV);
                nVox2Use_now = [];
            end

            %% run each classification scheme

            for gg=1:nGroups

                % taking out only trn and test trials with the relevant two
                % labels that we're trying to discriminate
                inds2use_tst = ismember(binLabs_tst, groups(gg,:));
                inds2use_trn = ismember(binLabs_trn, groups(gg,:));

                group_labs_tst = binLabs_tst(inds2use_tst);
                group_labs_trn = binLabs_trn(inds2use_trn);
                dat_this_group_tst = squeeze(allDat_tst(inds2use_tst, :));
                dat_this_group_trn = squeeze(allDat_trn(inds2use_trn, :));
                cv_this_group_tst = cvLabs_tst(inds2use_tst);
                cv_this_group_trn = cvLabs_trn(inds2use_trn);
                
                runs_using_trn = unique(cv_this_group_trn);
                assert(all(runs_using_trn==unique(cv_this_group_tst)))
                % train/test the decoder - using custom code to
                % cross-validate (same procedure used for within-cond
                % decoding)
                [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(dat_this_group_trn,group_labs_trn,...
                    cv_this_group_trn,dat_this_group_tst, group_labs_tst,...
                    cv_this_group_tst,class_str,100,nVox2Use_now,voxStatTable(:,runs_using_trn),0);

                % compute accuracy on the subset of test trials in each
                % condition separately                  
                allacc(vv,cc,gg) = mean(thesepredlabs==group_labs_tst);
                alld(vv,cc,gg) = get_dprime(thesepredlabs,group_labs_tst,unique(group_labs_tst));

                % confidence is the distance to incorrect - distance to
                % correct. want a positive number (far from incorrect)
                conf = normEucDist(:,2) - normEucDist(:,1);
                conf(group_labs_tst==groups(gg,2)) = -conf(group_labs_tst==groups(gg,2));
                % check these confidence labels to make sure they track -
                % always positive when classifier is correct, negative when
                % classifier makes a mistake.
                assert(all(conf(thesepredlabs==group_labs_tst)>0) && all(conf(thesepredlabs~=group_labs_tst)<0))
                conf_this_cond(inds2use_tst) = conf;

            end
           
            % put confidence labels back together over all conditions
            allconf(vv, condLabs==cc) = conf_this_cond;
        end
    end
    
    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','rt','correct','condlabs','allacc_rand','alld_rand');

end