function OriSpin2_unwarp()
% recon does unwarping of reconstructed multiband data.
% It looks for nifti in inpath/nifti and for DICOM in inpath/dicom (you
% can change these defaults if you want by changing the strings that are
% entered for nii_dir and dcm_dir line 39 and 40, respectively). The 
% function writes unwarped .nii.gz files to inpath/nifti. The moco flag is
% off by default, however, if you want the CFMRI routine to do moco (using
% AFNI) for you before the unwarping, you can opt into this by setting this
% flag to 1. Finally, this function also writes a log files to your 
% outpath, specifying the files that went into each unwarping call.
%
% inpath =
% '/mnt/neurocube/local/serenceslab/Rosanne/WM_Distract/DataRaw/S04/Session1'

inpath = '/usr/local/serenceslab2/maggie/OriSpin2/DataRaw/S07/Session4';
domoco = false;

%% catch errors
% If there is no outpath, write to the inpath
if nargin < 2, domoco = false; end
% If there are no file seperators ('/') at the end of your in or oupaths, put one there:
if ~strcmp(inpath(end), filesep)
    inpath = [inpath filesep];
end



%% Set defaults
my_path = pwd; % I'm here
% topup scripts path. Notice how I have my own copy of all the scripts...
% this is just a precaution. A tmp folder will be made wherever the script
% is ran from, and in this function I cd to the outpath, so a bunch of junk 
% outout should go there. But just in the unlikely case stuff does get 
% written to the path where the topup scrips live, two people doing recon/
% unwarping at the same time, might interfere in unexpected ways. 
topup_dir = '/mnt/neurocube/local/serenceslab/maggie/mFiles/run_unwarping/'; 
addpath(genpath(topup_dir));    % topup scripts path
nii_dir = 'Niftis';             % my reconstructed will be found under inpath/nii_dir
dcm_dir = 'Dicoms';             % my dicom will be found under inpath/dcm_dir
numcores = 8;
if domoco, outstr = '_topup_moco'; mocostr = ' -domoco'; else outstr = '_topup'; mocostr = ''; end


%% Read reconstructed Nifti folder names
niis = dir([inpath, nii_dir, '/*mbrecon.nii.gz']);
num_niis = size(niis,1);


%% Read sDir folder names (and run imseq if needed)
s_dirs = dir([inpath, dcm_dir, '/s*']);
num_dcms = NaN(1,length(s_dirs));
for s_idx = 1:length(s_dirs)
    if ~isempty(dir([inpath, dcm_dir, '/', s_dirs(s_idx).name, '/*MRDC*'])) %see if there are MRDC files
        fprintf(['\ndicom still need to be reordered. Running ''imseq'' now\n']);
        unix(['imseq ', inpath, dcm_dir, '/', s_dirs(s_idx).name]);
    end
    num_dcms(s_idx) = numel(dir([inpath, dcm_dir, '/', s_dirs(s_idx).name, '/*CFMRI*']));
end
% find topups
topup_dirs = s_dirs(num_dcms==144);
% find functionals
[~,sortingInds] = sort(num_dcms,'descend');
s_dirs = s_dirs(sort(sortingInds(1:num_niis)));


%% Check if unwarped files already exist before starting
do_topup = ones(num_niis,1);
for file_idx = 1:num_niis
    if exist([inpath, 'Niftis/', niis(file_idx).name(1:end-7), outstr, '.nii.gz'],'file')
        rsp =  input(['A topupped nifti (',niis(file_idx).name(1:end-7), outstr, '.nii.gz) already exists. Overwrite? (y/n) '], 's');
        if strcmp(rsp, 'n')
            do_topup(file_idx) = 0;
        else
            unix(['rm ', niis(file_idx).name(1:end-7), outstr, '.nii.gz']);
        end
    end
end
    
   

%% Do unwarping   
topuplogfile = [inpath, 'Niftis/topup.log'];
fid = fopen(topuplogfile,'w');
fprintf(fid,'Files used for topup\n.nii.gz\t\t\t\tfwd\trvs\n');
for file_idx = 1:num_niis
    if do_topup(file_idx)
        % log the input files to my recon call
        fprintf(fid,['\n', niis(file_idx).name, '\t', topup_dirs(1).name, '\t', topup_dirs(2).name]);
        % if they don't exist in the outpath yet, move my topup sDirs there
        cd([inpath, 'Niftis']) % <--important, needs all its shit in one place in order to work
        for n = 1:2
            if ~exist(topup_dirs(n).name, 'dir')
                unix(['cp -r ', inpath, dcm_dir, '/', topup_dirs(n).name, ' ', inpath, nii_dir],'-echo');
            end    
        end
        % and if a tmp folder is in my outpath, get rid of it
        if exist([inpath, nii_dir, 'tmp'],'dir')
            unix(['rm -rf ', inpath, nii_dir, 'tmp'],'-echo');
        end
        % actual unwarping
        file_to_unwarp = niis(file_idx).name;
        [s,~] = unix([topup_dir, 'run_topup_serences -d1 ', topup_dirs(1).name,' -d2 ', topup_dirs(2).name,...
            ' -i ', inpath, nii_dir, filesep, file_to_unwarp(1:end-7), ' -o ', inpath, nii_dir, filesep, file_to_unwarp(1:end-7), outstr, mocostr], '-echo');
        % check for error
        if s > 0 
            fprintf(['\nUnwarping (a.k.a. toptup) failed for ', file_to_unwarp, '.nii.gz\n\n'])
        else % if no error I'm gonna clean up
            stupid_logs = dir(['topup_*.log']);
            if ~isempty(stupid_logs)
                for numlogs = 1:length(stupid_logs)
                    unix(['rm ', stupid_logs(numlogs).name],'-echo');
                end
            end
        end
    end
end
fclose(fid); 
for n = 1:2 % and remove my sDirs again too
    if exist(topup_dirs(n).name, 'dir')
        unix(['rm -r ', inpath, nii_dir, filesep, topup_dirs(n).name],'-echo');
    end
end
   

%% return to where I started
cd(my_path)
   
