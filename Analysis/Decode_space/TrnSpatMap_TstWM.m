%% Spatial decoding analysis
% Train and test linear decoder, using data from independent spatial
% working memory mapping task as training set and testing within each
% condition of main task. 
% grouping spatial positions into 8 bins 45 deg wide, then doing binary
% decoding between bins that are 180 deg apart (e.g. 0 vs 180). 
% Saves the results in a mat file, can plot it using a separate script
% (plotTrnSpatMap_TstWM.m)
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
rndseed = 324454;
rng(rndseed,'twister');

% define the spatial position bins
nBins=8;
bin_centers=[0:45:359];
bin_size=diff(bin_centers(1:2));
% going to do 4 different classifications - bin 1 versus bin 5, 2
% versus 6, etc.
groups = [(1:4)',(5:8)'];
nGroups = size(groups,1);
nConds=2;

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('SpatLocSignalByTrial_%s.mat',substr));
    load(fn2load);
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnSpatMap_TestWM_BigBins_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
    %% loop over ROIs and run the model for each.

    for vv = 1:length(ROI_names) 
         %% check if there's data here

        if length(locSig)<vv || isempty(locSig(vv).dat_avg) || size(locSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s\n',ROI_names{vv});
        
        %% pull out the data for localizer (training), and main task (testing)
        trnPosLabs_orig = locSig(vv).StimPos;
        % IMPORTANT: SWITCH FROM THE COORDINATE SYSTEM OF SPATIAL LOCALIZER
        % TO THE COORDINATE SYSTEM OF MY SPATIAL WM TASK
        % localizer positions go clockwise from top/vertical, while spatial
        % WM positions go counter-clockwise from rightmost pt.
        trnPosLabs = mod(90 - trnPosLabs_orig,360);
        wedge_width_deg = 15;   % how many degrees wide is the wedge?
        left_wedge_edges = unique(trnPosLabs)-wedge_width_deg/2;
        right_wedge_edges = unique(trnPosLabs)+wedge_width_deg/2;
        
        trnDat = locSig(vv).dat_avg;
        % subtract mean over voxels
        trnDat = trnDat - repmat(mean(trnDat,2), 1, size(trnDat,2));

        tstPosLabs = mainSig(vv).targPos;
        tstDat = mainSig(vv).dat_avg;
        % subtract mean over voxels
        tstDat = tstDat - repmat(mean(tstDat,2), 1, size(tstDat,2));
         
        condLabs = mainSig(vv).condLabs;
        dist_to_bound = mainSig(vv).dist_to_real_bound;
        
        if vv==1
            % preallocate array here
            allacc = nan(length(ROI_names), nConds, nGroups);
            alld = nan(length(ROI_names), nConds, nGroups);
            allconf = nan(length(ROI_names), numel(tstPosLabs));
            allacc_rand = nan(length(ROI_names), nConds, nGroups, nPermIter);
            alld_rand = nan(length(ROI_names), nConds, nGroups, nPermIter);
        end
        
        % bin these for classifier - want 8 bins that are roughly centered at
        % 0, 45, 90, 135. 
        % For the training set here (spatial mapping task) there were 24
        % total "wedge" positions - here i'm including 4 wedge positions
        % per bin. So, 25% of the wedges which are on a border between
        % bins, end up getting used twice for the two diff decoders.
        
        trnBinLabs = zeros(size(trnPosLabs,1),2);   %make a second column because some wedges belong to two bins 
        for bb=1:nBins
            inds_this_bin = abs(trnPosLabs-(bin_centers(bb)))<wedge_width_deg*2 | abs((trnPosLabs-360)-(bin_centers(bb)))<wedge_width_deg*2;
            already_counted = trnBinLabs(:,1)~=0;
            trnBinLabs(inds_this_bin & ~already_counted,1) = bb;
            trnBinLabs(inds_this_bin,2) = bb;
        end
        
        % now going in and removing the trials with the in-between wedge
        % positions, because won't use those in decoding.
%         inds2use_trn = trnBinLabs~=0;
%         trnBinLabs = trnBinLabs(inds2use_trn);
%         trnDat = trnDat(inds2use_trn,:);
        
        assert(~any(trnBinLabs(:)==0))
%         neach_trn = sum(repmat(trnBinLabs,1,nBins)==repmat((1:nBins),size(trnBinLabs,1),1));
%         assert(all(neach_trn==neach_trn(1)));
        
        % same thing for test labels
        tstBinLabs = zeros(size(tstPosLabs));        
        for bb=1:nBins
            inds_this_bin = abs(tstPosLabs-(bin_centers(bb)-0.0001))<bin_size/2 | abs((tstPosLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
            tstBinLabs(inds_this_bin) = bb;
        end
        assert(~any(tstBinLabs==0))
        neach_tst = sum(repmat(tstBinLabs,1,nBins)==repmat((1:nBins),size(tstBinLabs,1),1));

        %% voxel selection from the training set 
       
        if ~isempty(nVox2Use) && nVox2Use<size(trnDat,2)
            fprintf('running voxel selection f-test for %s %s\n',substr, ROI_names{vv})
            
            pvals = zeros(size(trnDat,2), 1);
            dat = trnDat;
            lab = trnBinLabs(:,1);
            parfor vx = 1:size(dat,2)
                 % choose the voxels        
               [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
               pvals(vx) = pvalue;
            end           
            [~, ind] = sort(pvals, 'ascend'); 
            vox2use = ind(1:nVox2Use);
        else            
            vox2use = 1:size(trnDat,2);
        end
        
        %% run each classification scheme
        
        for gg=1:nGroups
            
            % taking out only trn and test trials with the relevant two
            % labels that we're trying to discriminate
            inds2use_trn = ismember(trnBinLabs(:,1), groups(gg,:)) | ismember(trnBinLabs(:,2), groups(gg,:));
            inds2use_tst = ismember(tstBinLabs, groups(gg,:));
            
            reallabstrn = trnBinLabs(inds2use_trn,:);
            % make sure we grab whichever label is of interest for this
            % decoder, not the other one.
            use_second_column = ~ismember(reallabstrn(:,1),groups(gg,:));
            reallabstrn(use_second_column,1) = reallabstrn(use_second_column,2);
            reallabstrn = reallabstrn(:,1);
            assert(sum(reallabstrn==groups(gg,1))==sum(reallabstrn==groups(gg,2)))
            
            reallabstst = tstBinLabs(inds2use_tst,1);
            
            assert(all(ismember(reallabstrn,groups(gg,:))));
            assert(all(ismember(reallabstst,groups(gg,:))));
            
            trnDatUse = trnDat(inds2use_trn, vox2use);
            tstDatUse = tstDat(inds2use_tst, vox2use);
            
            % train/test the decoder
            [thesepredlabs,normEucDist] = normEucDistClass(trnDatUse,tstDatUse,reallabstrn);
            
            % compute accuracy on the subset of test trials in each
            % condition separately
            
            tstconds = condLabs(inds2use_tst);            
            for cc=1:nConds
                allacc(vv,cc,gg) = mean(thesepredlabs(tstconds==cc)==reallabstst(tstconds==cc));
                alld(vv,cc,gg) = get_dprime(thesepredlabs(tstconds==cc),reallabstst(tstconds==cc),unique(reallabstst));               
            end
            
            % confidence is the distance to incorrect - distance to
            % correct. want a positive number (far from incorrect)
            conf = normEucDist(:,2) - normEucDist(:,1);
            conf(reallabstst==groups(gg,2)) = -conf(reallabstst==groups(gg,2));
            % check these confidence labels to make sure they track -
            % always positive when classifier is correct, negative when
            % classifier makes a mistake.
            assert(all(conf(thesepredlabs==reallabstst)>0) && all(conf(thesepredlabs~=reallabstst)<0))
            allconf(vv,inds2use_tst) = conf;
                  
            fprintf('%s %s gg=%d, performance on real data is %.2f and %.2f, starting random shuffles over %d iters...\n',...
                substr,ROI_names{vv},gg,allacc(vv,1,gg),allacc(vv,2,gg),nPermIter)
            
            % now doing the permutation test, shuffle labels 1000 times.
            randaccs_cond1 = nan(nPermIter, 1);
            randaccs_cond2 = nan(nPermIter, 1);
            randd_cond1 = nan(nPermIter, 1);
            randd_cond2 = nan(nPermIter, 1);
            
            % doing the shuffling before parfor loop 
            randlabs_all = zeros(size(reallabstrn,1),nPermIter);
            for ii=1:nPermIter
                randlabs_all(:,ii) = reallabstrn(randperm(numel(reallabstrn)));
            end          
            parfor ii=1:nPermIter
                % randomize training set labels
                randlabstrn = randlabs_all(:,ii);                
                % run classifier with the random labels
                [thesepredlabs,~] = normEucDistClass(trnDatUse,tstDatUse,randlabstrn);
                % get performance in each condition, for the random decoder
                randaccs_cond1(ii) = mean(thesepredlabs(tstconds==1)==reallabstst(tstconds==1));
                randaccs_cond2(ii) = mean(thesepredlabs(tstconds==2)==reallabstst(tstconds==2));
                randd_cond1(ii) = get_dprime(thesepredlabs(tstconds==1),reallabstst(tstconds==1),unique(reallabstst));
                randd_cond2(ii) = get_dprime(thesepredlabs(tstconds==2),reallabstst(tstconds==2),unique(reallabstst));
            end
           
            % put everything into a big array for saving
            allacc_rand(vv,1,gg,:) = randaccs_cond1;
            allacc_rand(vv,2,gg,:) = randaccs_cond2;
            alld_rand(vv,1,gg,:) = randd_cond1;
            alld_rand(vv,2,gg,:) = randd_cond2;
        end
        
    end
    
    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','rt','correct','condlabs','allacc_rand','alld_rand');

end