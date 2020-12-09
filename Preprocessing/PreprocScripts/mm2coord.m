function [inds] = mm2coord(mm,sform)
% convert mm (world coords) to coordinates in matrix (voxel coords)
% sform is a field in nifti header, nifti.sform

assert(size(mm,2)==3);

% this point defines the mm coordinates of the corner voxel, having index
% (1,1,1). this is the center of it.
corner_mm = repmat(sform(1:3,4)',size(mm,1),1);

% this matrix describes a mapping from mm to coords
map = pinv(sform(1:3,1:3));

inds = (mm - corner_mm) *map';

% add a unit to account for zero-indexing
inds = inds + 1;