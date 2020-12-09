function status = Label2VOI(sbj, ROIs, hemis, func_template, regfile, outdir, subjects_dir, ROINames)
% Converts all label files for ROIs <ROIs> and hemispheres <hemis> for
% subject <sbj> to 3D NIFTI binary masks, using <regfile> and
% <func_template> for alignment.
%
% The optional argument ROINames allows the user to specify output names
% that are different from the label names used as input. This may come in
% handy when the input name is overly long.

status = 0;

if nargin < 8, ROINames = ROIs; end
if nargin < 7, subjects_dir = []; end
if nargin < 6, outdir = './'; end
if nargin < 5
    out('\nERROR: Not enough input arguments.\n')    
    status = 1;
    return
end


if isempty(subjects_dir)
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
    if isempty(SUBJECTS_DIR)
        error('FreeSurfer subject directory undefined.');        
    else
        subjects_dir = SUBJECTS_DIR;
    end
end

for hemi_idx = 1:length(hemis)
    for roi_idx = 1:length(ROIs)
        labelfile = [subjects_dir, filesep, sbj, '/label/labels/', hemis{hemi_idx}, '.', ROIs{roi_idx}, '.label'];
        if exist(labelfile, 'file')
            unixcmd = ['mri_label2vol --label ', labelfile, ' --temp ', func_template, ' --subject ', sbj, ' --hemi ', hemis{hemi_idx}, ...
                ' --o ', outdir, hemis{hemi_idx}, '_', ROINames{roi_idx}, '.nii.gz --proj frac 0 1 0.1 --fillthresh 0.3 --reg ', regfile];
            [s,~] = unix(unixcmd, '-echo');
            if s==1
                out(['\nERROR: An error occured while processing label ', labelfile, '. Exiting.\n']);
                status = 1;
                return
            end                
        end            
    end %roi_idx   
end



end