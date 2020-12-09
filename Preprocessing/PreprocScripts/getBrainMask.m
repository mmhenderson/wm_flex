% Make a brain mask file for your subject in the functional space
% corresponding to this experiment. This script finds a brain mask file
% that lives in Doreti folder in anatomical space, then uses your
% regheadercenterMOD.mat transformation matrix to map it into the
% functional space that your preprocessed niftis are in. Can then use it to
% mask out good voxels for further analyses!

% the resulting file is saved in DataPreproc/S##/BrainMask_REG2FUNC.nii

% MMH 8/17/18

%%
clear
close all

sub = 'S03';
sub_doreti = 'CI';
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

doreti_path = '/mnt/neurocube/local/serenceslab/Doreti/';

unzip = 0;

%%

if ~exist([exp_path 'DataPreproc/' sub '/regheadercenterMOD.mat'],'file')
    error('need your registration matrix first!')
end

% Need a reference volume, this should be the first volume from your
% functional data (has to be session 1!!)
% reference_volume = [exp_path, 'DataPreproc/' sub '/MCTemplate01.nii.gz'];
reference_volume = [exp_path, 'DataPreproc/' sub '/MCTemplateXFM01.nii.gz'];

subpath = [exp_path, 'DataPreproc/', sub, '/']; % path to this subject's preprocessed data

% get the brain mask file into nii format
file2align = [doreti_path 'ANAT/'  sub_doreti '/mri/brainmask.nii.gz'];
if ~exist(file2align,'file')
    err = unix(['mri_convert ' doreti_path 'ANAT/'  sub_doreti '/mri/brainmask.mgz ' file2align]);
    if err
        error('your mri_convert command failed!')
    end
end

file2save = [exp_path 'DataPreproc/' sub '/BrainMask_REG2FUNC.nii.gz'];

fprintf('loading your original brain mask from %s\n',file2align)

fprintf('saving your transformed brain mask to %s\n',file2save)
%%
% make the transformation
err =  unix(['flirt -in ', file2align ' -ref ', reference_volume,...
    ' -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', file2save]);
if err
    error('your flirt command failed!')
else
    disp('SUCCESS')
end


% and unzip it if you want to
if unzip
    fprintf('unzipping your transformed brain mask...\n')

    err = unix(['gzip -dk ' file2save]);
    if err
        error('your unzip command failed!')
    else
        disp('SUCCESS')
    end
end