%% Spatial decoding analysis
% Train and test linear decoder, using data from independent spatial
% working memory mapping task as training set and testing within each
% condition of main task. 
% grouping spatial positions into 8 bins 45 deg wide, then doing binary
% decoding between bins that are 180 deg apart (e.g. 0 vs 180). 
% Saves the results in a mat file, can plot it using a separate script
% (plotTrnSWM_TstWM.m)
%%
clear
close all;

% sublist = [2:7];
sublist = [3];
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

    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnTestSpatMap_%s_max%dvox_%s.mat',class_str, nVox2Use,substr));
    
    %% loop over ROIs and run the model for each.

    for vv = 1:length(ROI_names) 
         %% check if there's data here

        if length(locSig)<vv || isempty(locSig(vv).dat_avg) || size(locSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s\n',ROI_names{vv});
        
        %% pull out the data for localizer (training), and main task (testing)
        trnPosLabs = locSig(vv).StimPos;
        wedge_width_deg = 15;   % how many degrees wide is the wedge?
        left_wedge_edges = unique(trnPosLabs)-wedge_width_deg/2;
        right_wedge_edges = unique(trnPosLabs)+wedge_width_deg/2;
        
        trnDat = locSig(vv).dat_avg;
        % subtract mean over voxels
        trnDat = trnDat - repmat(mean(trnDat,2), 1, size(trnDat,2));

        runLabs = locSig(vv).runLabs;

        if vv==1
            % preallocate array here
            allacc = nan(length(ROI_names), nGroups);
            alld = nan(length(ROI_names), nGroups);
            allconf = nan(length(ROI_names), numel(trnPosLabs));
            allacc_rand = nan(length(ROI_names), nGroups, nPermIter);
            alld_rand = nan(length(ROI_names), nGroups, nPermIter);
        end
        
        % bin these for classifier - want 8 bins that are roughly centered at
        % 0, 45, 90, 135. 
        % For the training set here (spatial mapping task) there were 24
        % total "wedge" positions - to have bins that are centered at the
        % right orients, we can ignore some of the bins. Basically taking
        % the two bins right and left of 0, the two bins right and left of
        % 45, etc. the one that falls in between (15-30 etc) gets ignored
        % here. alternative is to include that bin in both decoders but
        % this would lose spatial specificity.
        
        trnBinLabs = zeros(size(trnPosLabs));
        for bb=1:nBins
            inds_this_bin = abs(trnPosLabs-(bin_centers(bb)))<wedge_width_deg | abs((trnPosLabs-360)-(bin_centers(bb)))<wedge_width_deg;
            trnBinLabs(inds_this_bin) = bb;
        end
        
        % now going in and removing the trials with the in-between wedge
        % positions, because won't use those in decoding.
        trials2use_trn = trnBinLabs~=0;
        trnBinLabs = trnBinLabs(trials2use_trn);
        trnDat = trnDat(trials2use_trn,:);
        runLabs = runLabs(trials2use_trn);
        nCV = numel(unique(runLabs));
        
        assert(~any(trnBinLabs==0))
        neach_trn = sum(repmat(trnBinLabs,1,nBins)==repmat((1:nBins),size(trnBinLabs,1),1));
        assert(all(neach_trn==neach_trn(1)));
        
       
        %% voxel selection from the training set 
       
        if ~isempty(nVox2Use) && nVox2Use<size(trnDat,2)
            fprintf('running voxel selection f-test for %s %s\n',substr, ROI_names{vv})
            voxStatTable = zeros(size(trnDat,2),nCV);
            for rr = 1:nCV
                inds = runLabs~=rr;
                pvals = zeros(size(trnDat,2), 1);
                dat = trnDat(inds,:);
                lab = trnBinLabs(inds,:);
                parfor vx = 1:size(trnDat,2)
                     % choose the voxels        
                   [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                   pvals(vx) = pvalue;
                end 
                voxStatTable(:,rr) = pvals;
            end
            nVox2Use_now = nVox2Use;
        else            
            % put in a placeholder here because using all voxels
            voxStatTable = zeros(size(trnDat,2),nCV);
            nVox2Use_now = [];
        end
        
        %% run each classification scheme

        for gg=1:nGroups
            
             % labels that we're trying to discriminate
            inds2use_trn = ismember(trnBinLabs, groups(gg,:));
            reallabstrn = trnBinLabs(inds2use_trn);
            runlabstrn = runLabs(inds2use_trn);
            trnDatUse = trnDat(inds2use_trn, :);

            runs_using = unique(runlabstrn);

            % train/test the decoder - use a custom function to do
            % cross-validation here. 
            [~,~,thesepredlabs,normEucDist] = my_classifier_cross_wconf(trnDatUse,reallabstrn,...
                runlabstrn,trnDatUse, reallabstrn,...
                runlabstrn,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),0);
 
            allacc(vv,gg) = mean(thesepredlabs==reallabstrn);
            alld(vv,gg) = get_dprime(thesepredlabs,reallabstrn,unique(reallabstrn));               

            % confidence is the distance to incorrect - distance to
            % correct. want a positive number (far from incorrect)
            conf = normEucDist(:,2) - normEucDist(:,1);
            conf(reallabstrn==groups(gg,2)) = -conf(reallabstrn==groups(gg,2));
            % check these confidence labels to make sure they track -
            % always positive when classifier is correct, negative when
            % classifier makes a mistake.
            assert(all(conf(thesepredlabs==reallabstrn)>0) && all(conf(thesepredlabs~=reallabstrn)<0))
            allconf(vv,inds2use_trn) = conf;
                  
            fprintf('%s %s gg=%d, performance on real data is %.2f, starting random shuffles over %d iters...\n',...
                substr,ROI_names{vv},gg,allacc(vv,gg),nPermIter)
            
            % now doing the permutation test, shuffle labels 1000 times.
            randaccs = nan(nPermIter, 1);           
            randd = nan(nPermIter, 1);
                       
            % doing the shuffling before parfor loop 
            randlabs_all = zeros(size(reallabstrn,1),nPermIter);
            for ii=1:nPermIter
                for cv=1:numel(runs_using)
                    % shuffle the data from one cross-validation fold at a time, so we
                    % don't un-balance the training sets. 
                    inds = runlabstrn==runs_using(cv);
                    dat2shuff=reallabstrn(inds);
                    randlabs_all(inds,ii) = dat2shuff(randperm(numel(dat2shuff)));
                end       
                
            end          
            parfor ii=1:nPermIter
                % randomize training set labels
                randlabstrn = randlabs_all(:,ii);                
                % run classifier with the random labels
                [~,~,thesepredlabs,~] = my_classifier_cross_wconf(trnDatUse,randlabstrn,...
                    runlabstrn,trnDatUse, randlabstrn,...
                    runlabstrn,class_str,100,nVox2Use_now,voxStatTable(:,runs_using),0);

                % get performance in each condition, for the random decoder
                randaccs(ii) = mean(thesepredlabs==reallabstrn);                
                randd(ii) = get_dprime(thesepredlabs,reallabstrn,unique(reallabstrn));
            end
           
            % put everything into a big array for saving
            allacc_rand(vv,gg,:) = randaccs;           
            alld_rand(vv,gg,:) = randd;
            
        end
        
    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld','allconf','allacc_rand','alld_rand');

end