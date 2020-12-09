% Make ROIs for visual (retinotopically-defined) ROIs in each subject, in
% the space needed for oriSpin.
clear
close all

% set inputs
FS_sbj = 'CP';
subnum = '07';
silent = 0;

% set paths (end all with filesep!)
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

feat_path = [exp_path, 'AnalyzeLocalizer/S', char(subnum), '/feats/AllSessionsFE.gfeat/cope1.feat/'];
out_path = [exp_path, 'Samples/'];
beh_path = [exp_path, 'DataBehavior/S', char(subnum), '/'];
func_path = [exp_path, 'DataPreproc/S', char(subnum), '/'];
hemis = {'lh', 'rh'};

%--> FS_sbj: needs to be the subject ID in Doreti.
%--> featdir: Directory where the FEAT output is (this is a .gfeat if
%    multiple localizer runs were done in a session, or simply a .feat dir
%    if only one localizer was ran 


%% Registration matrix
%The registration file is needed to know how voxels from your functional
%runs map onto voxels from the Doreti retinotopy scan. Thus: a functional
%template from your experiment needs to be registered to the orig.mgz. If
%no reg file exists, user must make one now or there's no point continuing.
fprintf('\n\n##### ONE MORE REGISTRATION #####\n\n')

%%% start with manual registration (i.e. tkregister2)
ManReg = true;
if exist([func_path, 'Func2Anat.dat'], 'file')
    rsp = input('A manual registration file already exists for this subject and session, do you want to overwrite? (y/n) ', 's');
    if strcmp(rsp,'n')
        ManReg = false;
    end
end
if ManReg
    [s,w] = unix(['tkregister2 --mov ' func_path 'MCTemplateXFM01.nii.gz --s ' char(FS_sbj) ' --regheader --reg ' func_path 'Func2Anat.dat']);
    if s == 1 % non-zero value returned by unix means the command has failed
        fprintf(['\n\nERROR: Failed to start manual registration. \n\n' 'ERROR MESSAGE FROM UNIX WRAPPER:' w '\n\n']);
        return
    end
end
%%% finetune manual with automatic registration
AutoReg = true;
if exist([func_path, 'Func2Anat_auto.dat'], 'file')
    rsp = input('An automatic registration file already exists for this subject and session, do you want to overwrite? (y/n) ', 's');
    if strcmp(rsp,'n')
        AutoReg = false;
    end
end
if AutoReg
    fprintf('\n...Performing automatic registration using boundary based alignment, give it a minute...\n');
    unix(['bbregister --s ', char(FS_sbj), ' --mov ' func_path 'MCTemplateXFM01.nii.gz --reg ' func_path 'Func2Anat_auto.dat --init-reg ' func_path 'Func2Anat.dat --bold']);
    fprintf('Finished automatic registration!\n');
end
%%% Visual check of registration
check_auto = input('Do you want to visually inspect the result of the automatic registration? (y/n) ', 's');
if strcmp(check_auto,'y')
    unix(['tkregister2 --mov ' func_path 'MCTemplateXFM01.nii.gz --reg ' func_path 'Func2Anat_auto.dat --surf']);
end


%% Create Masks
% Check for the existence of a VOI folder. If it doesn't exist, that means
% we need to create one and convert labels to masks. If it does exist, 
% prompt for overwrite. Labels we have from the retinotopy, they're labels
% on the surface representing different visual areas. We want to make these
% into volumes so that we know which voxels in our 4D functional data 
% correspond to which visual area. 

fprintf('\n\n##### CREATE RETINOTOPIC MASKS (AKA MAKE VOIs) #####\n\n')

% check to see if masks are already made or not
VisualMask = true;
if exist([exp_path, 'VOIs/S', char(subnum)], 'dir')
    if ~isempty(dir([exp_path, 'VOIs/S', char(subnum), '/*h_*nii.gz'])) % if there are visual masks found
        rsp = input('Visual area masks already exist in this subject''s VOI folder. Overwrite? (y/n) ', 's');
        if strcmp(rsp, 'n')
            VisualMask = false;
        end
    end
else
    mkdir([exp_path, 'VOIs/S', char(subnum)]);
end

% make Visual Area masks
if VisualMask
%     delete([experiment_path, 'VOIs/S' char(subnum) '/*h_*nii.gz']); %delete retino content from VOI folder 
    ROIs = {'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v', 'V3AB', 'hV4', 'LO1', 'LO2', 'IPS0', 'IPS1', 'IPS2', 'IPS3'};
    % Now convert the labels to volumes (or VOI's) with homespun function
    % LABEL2VOI(sbj, ROIs, hemis, func_template, regfile, outdir)
    s = Label2VOI(FS_sbj, ROIs, hemis, [func_path, 'MCTemplateXFM01.nii.gz'],...
        [func_path, 'Func2Anat_auto.dat'], [exp_path, 'VOIs/S' char(subnum) '/']);
    if s == 1
        out('\nERROR: Label conversion failed. Exiting.\n');
        return
    end
end

