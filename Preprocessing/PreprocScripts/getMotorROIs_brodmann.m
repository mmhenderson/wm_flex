% use the activation from digit localizer task to define motor-related ROIs
% in each subject.

% using the definition of brodmann areas done by freesurfer's recon-all
% function to define the rough boundaries of where areas should be, then these
% get thresholded with the digit localizer and saved as final ROIs.

% ROIs will be thresholded using contrast for the contralateral
% hemi (e.g. left M1 is defined using right>left contrast)

% MMH 11/17/20
%%
clear
close all

subnum = '07';
subinit = 'CP';

% set this to 1 if you want to make all the volumes in functional space.
% takes a while and only has to be done once per subject.
doVolumes = 1;

% these two contrasts have to correspond to copes 1 and 2!
contrasts = {'left>right','right>left'};
hemis = {'lh','rh'};
nHemis = length(hemis);
% set paths (end all with filesep!)

% defining the brodman areas that I want: areas 1-3 are primary
% somatosensory, area 4 is primary motor, ares 6 is premotor.
ba_parcels_to_use = {{'BA1','BA2','BA3a','BA3b'},{'BA4a','BA4p'},{'BA6'}};
final_roi_names = {'S1','M1','Premotor'};

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

retfolder = '/usr/local/serenceslab/Doreti/';
subretfolder = [retfolder 'ANAT/' subinit '/'];
func_path = [exp_path, 'DataPreproc/S', char(subnum), '/'];
% this extension gets appended to the final nifti files for each volume,
% indicating they were thresholded with digit localizer.
ROI_ext  = '_DIGITLOC_BA';

% this first dir is a place to store the un-thresholded VOIs after they are
% made, the second dir is our usual VOI folder
VOIdir1 = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
VOIdir2 = [exp_path, 'VOIs/S' char(subnum) '/'];


minVox = 1;

currentdir = pwd;

% load the annotation file so that we can get a list of all the names (this
% can really be any subject or hemisphere because they're all the same
% names)
% currently the way this is done is to make a volume for every single
% cortical area, then see how many intersect with the digit localizer.
annotfile = fullfile(retfolder, 'ANAT',subinit,'label','lh.BA_exvivo.annot');
[inds,labs,info]= read_annotation(annotfile);
allnames = info.struct_names(2:end);
allparcelnames = allnames;
parcels_4unix = allparcelnames;

for rr=1:length(allparcelnames)
    if contains(allparcelnames{rr},'&')
        % put an extra slash in here so that unix won't be confused by the & signs
        newname = allparcelnames{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '\' newname(andind:end)];
        parcels_4unix{rr} = newname;
        newname = allparcelnames{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '_and_' newname(andind+1:end)];
    
    end
end
clear info



%% first go into the Doreti folder, and make labels based on the anatomy

if doVolumes
    
    % check to see if labels are already made or not
    makelabs = true;
    if exist([subretfolder, 'label/parcels_all'], 'dir')
        if ~isempty(dir([subretfolder, 'label/parcels_all/*.label'])) % if there are labels in here already
            rsp = input('Labels already exist in this subject''s parcels_all folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makelabs = false;
            end
        end
    else
        mkdir([subretfolder, 'label/parcels_all']);
    end
    
    if makelabs
        cd(retfolder)
        for hh=1:nHemis
            unix(['mri_annotation2label --annotation BA_exvivo --subject ' subinit ' --hemi ' hemis{hh} ' --outdir ' subretfolder 'label/parcels_all/']);
        end
        cd(currentdir)
    end
    %% now turn them into VOI masks in functional space
    
    % here we are making a new folder to hold ROIs that were generated using
    % fsl parcellation scheme. These have to be processed a bit more before
    % they can be put into our normal pipeline.
    
    % check to see if masks are already made or not
    makemasks = true;
    if exist([exp_path, 'VOIs_parcel/S', char(subnum)], 'dir')
        if ~isempty(dir([exp_path, 'VOIs_parcel/S', char(subnum), '/*h_*nii.gz'])) % if there are visual masks found
            rsp = input('Visual area masks already exist in this subject''s VOIs_parcel folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makemasks = false;
            end
        end
    else
        mkdir([exp_path, 'VOIs_parcel/S', char(subnum)]);
    end
    
    % make Visual Area masks
    if makemasks
        delete([exp_path, 'VOIs_parcel/S' char(subnum) '/*h_*nii.gz']); %delete retino content from VOI folder
        %     ROIs_orig = {'S_precentral-sup-part','S_precentral-inf-part','S_central','S_postcentral',...
        %         'G_precentral','G_postcentral','G_parietal_sup',...
        %         'G&S_subcentral','G_pariet_inf-Supramar','G&S_paracentral','G_front_sup'};
        
        func_template = [func_path, 'MCTemplateXFM01.nii.gz'];
        regfile = [func_path, 'Func2Anat_auto.dat'];
        outdir = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
        if ~isdir(outdir)
            mkdir(outdir)
        end
        if ~exist(regfile,'file')
            error('registration to anat not made yet! you need to run getVisualROIs.m first')
        end
        % below is in replacement of label2VOI function
        for hemi_idx = 1:length(hemis)
            for roi_idx = 1:length(parcels_4unix)
                labelfile = [subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', parcels_4unix{roi_idx}, '.label'];
                
                
                if exist([subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', allparcelnames{roi_idx}, '.label'], 'file')
                    
                    unixcmd = ['mri_label2vol --label ', labelfile, ' --temp ', func_template, ' --subject ', subinit, ' --hemi ', hemis{hemi_idx}, ...
                        ' --o ', outdir, hemis{hemi_idx}, '_', parcels_4unix{roi_idx}, '.nii.gz --proj frac 0 1 0.1 --fillthresh 0.3 --reg ', regfile];
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
    
    
end
%% now load the localizer, and mask out the parts we need
%
for cc=1:length(contrasts)
    
    % find stats map
    stats_path = [exp_path 'AnalyzeDigitLocalizer/S' char(subnum) '/feats/AllSessionsFE.gfeat/cope' num2str(cc) '.feat/stats/'];
    cd(stats_path)
    
    %% load my t-stat map for all voxels
    nifti = load_nifti('tstat1.nii.gz');
    vmpVals = nifti.vol; %Store t-values for each voxel index that is in my visual areas in 'vmpVals'
    
    clear nifti
    %% find my threshold
    % this is all copied from Find_t_thresh.m
    
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
    
    %% define my final mask
    
    vals2use = vmpVals>my_t;
    
    %% loop through cortical ROIs in the contralateral hemi only
    
    hh=3-cc;
    % first looping over the areas we want to create
    for pp=1:length(ba_parcels_to_use)
        
        % then within each area, loop over sub-regions that will be in it
        for vv=1:length(ba_parcels_to_use{pp})
            
            roi_name_load = ba_parcels_to_use{pp}{vv};
            
            fn2use = dir([VOIdir1, hemis{hh}, '*', roi_name_load, '*.nii.gz']);
            if numel(fn2use)==1
                                               
                thisVOI = load_nifti(fullfile(fn2use(1).folder, fn2use(1).name));
                
                if vv==1
                    fullVol = zeros(size(thisVOI.vol));
                end
                % taking any voxels that belong in any component of the
                % area (combining areas)
                fullVol = fullVol | thisVOI.vol;
                
            else
                error('area %s not found\n',allparcelnames{vv})
            end
            
        end
        
        % now intersect it with the localizer map
        maskedVol = fullVol & vals2use;
        %         fprintf('%s %s %s: found %d voxels for %s\n',subinit,hemis{hh}, ROIs{vv},sum(maskedVol(:)),contrasts{cc});
        
        if sum(maskedVol(:))>minVox
            fprintf('%s %s %s: found %d voxels for %s\n',subinit,hemis{hh}, final_roi_names{pp},sum(maskedVol(:)),contrasts{cc});
            newVOI = thisVOI;
            newVOI.vol = maskedVol;
            savename = [VOIdir2 hemis{hh} '_' final_roi_names{pp} ROI_ext '.nii.gz'];
            err = save_nifti(newVOI,savename);
            if err~=0
                error('error saving new VOI file')
            end
        end
        
        
    end
    
    
    
end