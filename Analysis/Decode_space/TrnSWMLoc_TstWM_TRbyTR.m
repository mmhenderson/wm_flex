% MMH 10/29/20

% training/testing on SPATIAL WORKING MEMORY localizer, testing on task data.

clear
close all;

sublist = [7];
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
rndseed = 121323;
rng(rndseed,'twister');

nBins=8;
bin_centers=[0:45:359];
bin_size=diff(bin_centers(1:2));
% going to do 4 different classifications - bin 1 versus bin 5, 2
% versus 6, etc.
groups = [(1:4)',(5:8)'];
nGroups = size(groups,1);
nConds=2;

for ss=1:length(sublist)
    
    allchanresp = [];

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('SWMLocSignalByTrial_%s.mat',substr));
    load(fn2load);
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_TRbyTR_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
    %% loop over ROIs and run the model for each.

    for vv = 1:length(ROI_names) 
         %% check if there's data here

        if length(locSig)<vv || isempty(locSig(vv).dat_avg) || size(locSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s\n',ROI_names{vv});
        
        %% pull out the data for localizer (training), and main task (testing)
        trnPosLabs = locSig(vv).TargPos;
        trnDat = locSig(vv).dat_avg;
        % trying to get rid of baseline shifts here.
        trnDat = trnDat - repmat(mean(trnDat,2), 1, size(trnDat,2));

        tstPosLabs = mainSig(vv).targPos;
        tstDat = mainSig(vv).dat_by_TR; % [nTrials x nTRs x nVox]
        % trying to get rid of baseline shifts here. do subtraction within
        % a TR only. 
        tstDat = tstDat - repmat(mean(tstDat,3), 1, 1, size(tstDat,3));
%         
        condLabs = mainSig(vv).condLabs;

        if vv==1
            % preallocate array here
            nTRs=size(tstDat,2);
            allacc = nan(length(ROI_names), nConds, nGroups, nTRs);
            alld = nan(length(ROI_names), nConds, nGroups, nTRs);
            allconf = nan(length(ROI_names), numel(tstPosLabs), nTRs);
            allacc_rand = nan(length(ROI_names), nConds, nGroups, nTRs, nPermIter);
            alld_rand = nan(length(ROI_names), nConds, nGroups, nTRs, nPermIter);
        end
        
        % bin these for classifier - want 8 bins that are roughly centered at
        % 0, 45, 90, 135,...
        trnBinLabs = zeros(size(trnPosLabs));
        for bb=1:nBins
            inds_this_bin = abs(trnPosLabs-(bin_centers(bb)-0.0001))<bin_size/2 | abs((trnPosLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
            trnBinLabs(inds_this_bin) = bb;
        end
        assert(~any(trnBinLabs==0))
        neach_trn = sum(repmat(trnBinLabs,1,nBins)==repmat((1:nBins),size(trnBinLabs,1),1));
        assert(all(neach_trn==neach_trn(1)));
        
        % same thing for test labels
        tstBinLabs = zeros(size(tstPosLabs));        
        for bb=1:nBins
            inds_this_bin = abs(tstPosLabs-(bin_centers(bb)-0.0001))<bin_size/2 | abs((tstPosLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
            tstBinLabs(inds_this_bin) = bb;
        end
        assert(~any(tstBinLabs==0))
        
        
        %% voxel selection from the training set 
       
        if ~isempty(nVox2Use) && nVox2Use<size(trnDat,2)
            fprintf('running voxel selection f-test for %s %s\n',substr, ROI_names{vv})
            
            pvals = zeros(size(trnDat,2), 1);
            dat = trnDat;
            lab = trnBinLabs;
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
        
        for tr=1:nTRs
            
            for gg=1:nGroups

                % taking out only trn and test trials with the relevant two
                % labels that we're trying to discriminate
                inds2use_trn = ismember(trnBinLabs, groups(gg,:));
                inds2use_tst = ismember(tstBinLabs, groups(gg,:));

                reallabstrn = trnBinLabs(inds2use_trn);
                trnDatUse = squeeze(trnDat(inds2use_trn, vox2use));
                tstDatUse = squeeze(tstDat(inds2use_tst, tr,vox2use));
                
                % train/test the decoder
                [thesepredlabs,normEucDist] = normEucDistClass(trnDatUse,tstDatUse,reallabstrn);

                % compute accuracy on the subset of test trials in each
                % condition separately
                reallabstst = tstBinLabs(inds2use_tst);
                tstconds = condLabs(inds2use_tst);            
                for cc=1:nConds
                    allacc(vv,cc,gg,tr) = mean(thesepredlabs(tstconds==cc)==reallabstst(tstconds==cc));
                    alld(vv,cc,gg,tr) = get_dprime(thesepredlabs(tstconds==cc),reallabstst(tstconds==cc),unique(reallabstst));
                end
                               
                % confidence is the distance to incorrect - distance to
                % correct. want a positive number (far from incorrect)
                conf = normEucDist(:,2) - normEucDist(:,1);
                conf(reallabstst==groups(gg,2)) = -conf(reallabstst==groups(gg,2));
                % check these confidence labels to make sure they track -
                % always positive when classifier is correct, negative when
                % classifier makes a mistake.
                assert(all(conf(thesepredlabs==reallabstst)>0) && all(conf(thesepredlabs~=reallabstst)<0))
                allconf(vv,inds2use_tst,tr) = conf;
                      
                fprintf('%s %s gg=%d, tr=%d, performance on real data is %.2f and %.2f, starting random shuffles over %d iters...\n',...
                    substr,ROI_names{vv},gg,tr,allacc(vv,1,gg,tr),allacc(vv,2,gg,tr),nPermIter)

                % now doing the permutation test, shuffle labels 1000 times.
                randaccs_cond1 = nan(nPermIter, 1);
                randaccs_cond2 = nan(nPermIter, 1);
                randd_cond1 = nan(nPermIter, 1);
                randd_cond2 = nan(nPermIter, 1);
                parfor ii=1:nPermIter
                    % randomize training set labels
                    randlabstrn = reallabstrn(randperm(numel(reallabstrn)));               
                    % run classifier with the random labels
                    [thesepredlabs,~] = normEucDistClass(trnDatUse,tstDatUse,randlabstrn);
                    % get performance in each condition, for the random decoder
                    randaccs_cond1(ii) = mean(thesepredlabs(tstconds==1)==reallabstst(tstconds==1));
                    randaccs_cond2(ii) = mean(thesepredlabs(tstconds==2)==reallabstst(tstconds==2));
                    randd_cond1(ii) = get_dprime(thesepredlabs(tstconds==1),reallabstst(tstconds==1),unique(reallabstst));
                    randd_cond2(ii) = get_dprime(thesepredlabs(tstconds==2),reallabstst(tstconds==2),unique(reallabstst));
                end

                % put everything into a big array for saving
                allacc_rand(vv,1,gg,tr,:) = randaccs_cond1;
                allacc_rand(vv,2,gg,tr,:) = randaccs_cond2;
                alld_rand(vv,1,gg,tr,:) = randd_cond1;
                alld_rand(vv,2,gg,tr,:) = randd_cond2;

            end
        end
      
    end

    rt=mainSig(1).RTLabs;
    correct = mainSig(1).RespActual==mainSig(1).CorrectResp;
    condlabs=mainSig(1).condLabs;
    
    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','rt','correct','condlabs','allacc_rand','alld_rand');

end