% Get rough anatomical definitions of MD network ROIs.
%
% load the anatomical definitions of a parcellation of the MD network 
% from Fedorenko, Duncan , & Kanwisher, 2013
% downloaded from imaging.mrc-cbu.cam.ac.uk/imaging/MDsystem. 
% This script will convert that nifti file into a new nifti in
% the same coordinate system of the MNI-average brain from freesurfer. 
% (all happens in TAL space). This includes UPSAMPLING the resolution of
% the volume from 2mm vox to 1mm vox. 
% Next step is to project this file into subject (functional) space using
% alignment between MNI-average brain and subject anatomical.
% (getMDROIS.m does this)
%
% only need to run this once, it'll put a file in the
% main preprocessing folder which can be used for all subs. 

% NOTE these ROIs didn't get used in any of our main analyses for this
% experiment.

% MMH 11/22/19
%% 
clear

root = pwd;
filesepinds = find(root==filesep);
exp_path = root(1:filesepinds(end-1));
curr_subfolder = root(filesepinds(end-1)+1:end);

% this is the name of the new file we will make here
parcels_xfm_nii = fullfile(exp_path, curr_subfolder, 'MDROI_XFM_TAL.nii.gz');

%% load the nifti with the MD ROIs masked out
% need to have a copy of this file in our current folder, this is the thing
% we downloaded from Fedorenko paper
parcels = load_nifti(fullfile(exp_path, curr_subfolder, 'MDROI.nii.gz'));
if isempty(parcels)
    error('need to have the file MDROI.nii.gz in your current directory');
end

% sform tells us a mapping from i,j,k, in the volume to x,y,z in mm
parcels_sform = parcels.sform;
parcel_dims = size(parcels.vol);

% size of the voxels (isotropic) in the original nifti file
voxmm_orig = parcels.pixdim(2);

% different numbers for each region in the MD network
labels = unique(parcels.vol(:));

%% load the MNI average brain
% this file is same for all subjects, it's a reference brain
tal_brain_orig_mgz = '/mnt/neurocube/local/freesurfer6/average/mni305.cor.mgz';
% make sure we have a copy of it in this folder, in nii format
tal_brain_orig_nii = fullfile(exp_path,curr_subfolder,'mni305.cor.nii.gz');
if ~exist(tal_brain_orig_nii,'file')
    err = unix(['mri_convert ' tal_brain_orig_mgz ' ' tal_brain_orig_nii]);
    if err
        error('your mri_convert command failed!')
    end
end
avg_brain = load_nifti(tal_brain_orig_nii);
avg_brain_dims = size(avg_brain.vol);
% sform tells us a mapping from i,j,k, in the volume to x,y,z in mm
avg_brain_sform = avg_brain.sform;

%% put all the parcels into the space of the MNI brain
new_masked_vol = zeros(avg_brain_dims);

for ll = 1:length(labels)
   
    label_inds = find(parcels.vol==labels(ll));
    [i,j,k]  = ind2sub(parcel_dims,label_inds);
    label_coords = [i,j,k];
    
    % converting from indices into the original matrix, into coordinates
    % in a standard mm (world-centered) coordinate system
    % need to use the first sform matrix to do this
    label_mm = coord2mm([i,j,k],parcels_sform);
    
    % now converting from mm coordinates in standard space, to voxel
    % coordinates in the NEW matrix (need to use the other sform matrix)
    label_coords_new = round(mm2coord(label_mm, avg_brain_sform));
        
    % need to up-sample the labels here - loop over each voxel in the
    % original (lower-res) brain, and put the label everywhere within a
    % small range of that voxel. This expands the labels slightly.
    for mm = 1:length(label_coords_new)
        
        % find the voxels within a small range of original center voxel
        % should be 27 things here if vox size is 2 mm originally
        xc = label_coords_new(mm,1)-voxmm_orig/2:label_coords_new(mm,1)+voxmm_orig/2;
        yc = label_coords_new(mm,2)-voxmm_orig/2:label_coords_new(mm,2)+voxmm_orig/2;
        zc = label_coords_new(mm,3)-voxmm_orig/2:label_coords_new(mm,3)+voxmm_orig/2;
        
        [xc,yc,zc] = meshgrid(xc,yc,zc);
        
        coords = [xc(:),yc(:),zc(:)];
        
        % make sure no indices are out of range
        out_of_bounds = any(coords>avg_brain_dims,2) |...
            any(coords<1,2);
        coords(out_of_bounds,:) = [];

        % then back to linear indices 
        inds = sub2ind(avg_brain_dims, xc(:),yc(:),zc(:));

        % label all of those voxels with the same label they had before
        new_masked_vol(inds) = labels(ll);

        
    end

    % double checking that this worked and we labeled all the voxels that
    % we needed to (plus more)
    label_inds_new = sub2ind(avg_brain_dims,label_coords_new(:,1),label_coords_new(:,2),label_coords_new(:,3));
    assert(all(new_masked_vol(label_inds_new)==labels(ll)));

end

%% save a new nifti with the MD definitions
new_nii = avg_brain;
new_nii.vol = new_masked_vol;
save_nifti(new_nii, parcels_xfm_nii);
fprintf('saved file to %s\n',parcels_xfm_nii);