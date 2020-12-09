%% Merge all spatial localizer maps (for the 24 unique positions) into one thresholded map.
% During FEAT analysis, I ran a GLM with 24 predictors, and 29 total
% contrasts to be evaluated. Contrasts 1 looks for voxels that are 
% responsive for all stimulus positions (turns up very few voxels) 
% contrasts 2-5 look for voxels responsive for one of the 4 quadrants. 
% Contrasts 6-29 are each testing for voxels responsive for one of the 24 
% stimulus positions. 
% To get all voxels that are selective for any of the stimulus positions,
% we will threshold each map and merge the thresholded "masks" into a 
% single mask that is saved out.
%%
clear 
close all

sublist = [7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

for ss = sublist
    
    substr = sprintf('S%02d',ss);
    
    % find stats map
    feat_path = fullfile(exp_path,'AnalyzeSpatLocalizer',substr,'feats','AllSessionsFE.gfeat');
    
    template_nifti = load_nifti(fullfile(feat_path,'cope1.feat','stats','tstat1.nii.gz'));
    
    merged_mask = zeros(size(template_nifti.vol));
    
    % merging these 24 contrasts corresponding to the 24 positions.
    contrasts2merge = [6:29];
    for pp = contrasts2merge
    
        % go into the dir that has the t-map for this predictor of
        % interest.
        stats_path = fullfile(feat_path,sprintf('cope%d.feat',pp), 'stats');
        fprintf('doing thresholding on %s\n',stats_path);
        cd(stats_path);
        % now threshold this map
         % get FE dof
        [~, my_FE_dof] = unix('fslstats tdof_t1 -M'); % gives the dof I need for fixed effects (which is the test I ran), output is a string

        % make log p-stats map
        unix(['ttologp -logpout logp1 varcope1 cope1 ', num2str(str2num(my_FE_dof))]);

        % convert from log to p-value map
        unix(['fslmaths logp1 -exp p1']);

        % do FDR on p-value map and get probability threshold
        [~,prob_thresh] = unix('fdr -i p1 -q 0.05');

        % go from my p-threshold back into my t-stat
        my_t = abs(tinv(str2num(prob_thresh(28:end)),str2num(my_FE_dof)));
    
        % now, load the t-map and apply the threshold.
        t_nifti = load_nifti(fullfile(stats_path,'tstat1.nii.gz'));
        
        above_thresh = t_nifti.vol>my_t;
        
        % add these values to the merged mask.
        merged_mask = merged_mask | above_thresh;
        
    end
    
    template_nifti.vol = merged_mask;
    savename = fullfile(feat_path,'AllPredsMergeThresh.nii.gz');
    fprintf('Saving result to %s\n',savename);
    save_nifti(template_nifti,savename);
    
end



