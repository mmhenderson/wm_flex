function [mm] = coord2mm(inds,sform)
% convert coords in matrix (voxel coords) to mm (world coords)
% sform is a field in nifti header, nifti.sform

% subtract a unit to account for zero-indexing
inds = inds-1;

assert(size(inds,2)==3);

% this point defines the mm coordinates of the corner voxel, having index
% (1,1,1)
corner_mm = repmat(sform(1:3,4)',size(inds,1),1);

% this matrix describes a mapping from coords to mm. 
map = sform(1:3,1:3);

mm = inds*map'+corner_mm;