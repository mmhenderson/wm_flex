%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0) Prepare everything I need to run preprocessing
% (1) Coregister each session to Doreti anatomy
% (2) Apply transformation on all functional data
% (3) Motion Correction of all functional data to a single template
% (4) Detrending all functional data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
inID = 'topup.'; % will use this string to look for raw functional data (nifti)

% what do you wanna do?
WannaCoregister = 1;
WannaTransform = 1;
WannaMC = 1;
WannaDetrend = 1;

sub_list_big = {'S07'};
sub_doreti_list_big = {'CP'};
all_sessions_big = {[1,2,3,4]};   % list all sessions that exist for this sub (get merged to runs.list together)
preproc_sessions_big = {[4]}; % list the ones you want to do right now (those not done yet)

for xx = 1:length(sub_list_big)

    sub = sub_list_big{xx};
    sub_doreti = sub_doreti_list_big{xx};
    all_sessions = all_sessions_big{xx};
    preproc_sessions = preproc_sessions_big{xx};

    %% step 0: put my ducks in a row
    subpath = [exp_path, 'DataPreproc/', sub, '/']; % path to this subject's preprocessed data
    if ~exist(subpath,'dir'); mkdir(subpath); end 

    % if we're remaking any of the runs.lists now, we need to also remake the
    % merged list. Otherwise we can skip that.
    make_list_big = 0;


    raw_path = cell(length(all_sessions),1);
    invols = cell(length(all_sessions),1);
    runnrs = cell(length(all_sessions),1);

    for session = all_sessions

        % set path to raw functional data (make sure it ends with a filesep)
        raw_path{session} = [exp_path, 'DataRaw/', sub, '/Session'  num2str(session), '/Niftis/']; %path to raw data
    %     raw_path{session} = [exp_path 'DataRaw/' sub '/Niftis/'];
        % find volumes to be preprocessed (at raw_path)
        invols{session} = dir([char(raw_path(session)), '*', inID, '*nii.gz']); %lists the files in 'raw_path' directory (can use wildcards) and saves them into struct 'invols'
        if isempty(invols{session})
            fprintf(['\n\nERROR: No compressed NIFTI volumes found at ', char(raw_path(session)), ' with ID string ', inID, '. returning.\n\n']);
            return;
        end

        if exist([char(raw_path(session)), 'runs.list'], 'file')
    %         make_list = 0;
            fprintf('A runs.list file already exists for Session %d\n',session)
            make_list = input('Want to over-write it? (y/n) \n--> ', 's');
            if strcmp(make_list,'y')
                unix(['rm ' char(raw_path(session)) 'runs.list'])
                make_list = 1;
                make_list_big = 1;
            elseif strcmp(make_list, 'n')
                make_list = 0;
            else
                error('need to enter either y or n when prompted')
            end
        else
            make_list = 1;
            make_list_big = 1;
        end

        % Make a 'runs.list' in the raw data directory for this session (if still necessary)
        if make_list
            fprintf('To make your runs.list file for session %d, indicate volumes you want preprocessed\n',session );
            %Create and open file 'runs.list'
            cd(char(raw_path(session)));
            runslist = fopen('runs.list','w');
            %Ask
            exp = input(' which functional volumes were main WM task? \n please indicate as follows: [1 2 3 4] etc \n --> ');
            swm = input(' which functional volumes were SWM task? \n please indicate as follows [1 2 3 4] etc \n --> ');
            dwm = input(' which functional volumes were DWM task? \n please indicate as follows [1 2 3 4] etc \n --> ');
            
            spatloc = input(' which functional volumes were spatial mapping localizer? \n please indicate as follows: [1 2 3 4] etc \n --> ');
            digloc = input(' which functional volumes were digit localizer (no delay)? \n please indicate as follows: [1 2 3 4] etc \n --> ');
            mdloc = input(' which functional volumes were MD localizer? \n please indicate as follows: [1 2 3 4]  etc\n --> ');
            all_vols = sort([exp,swm,dwm,spatloc,digloc,mdloc]);
            for vols = 1:length(all_vols)
                %which volume is this
                volume = sprintf('%03d',vols);
                %is this an exp of a loc volume?
                if sum(exp==all_vols(vols))==1; text = [volume ' stim'];
                elseif sum(swm==all_vols(vols))==1; text = [volume ' swm'];
                elseif sum(dwm==all_vols(vols))==1; text = [volume ' dwm'];
                elseif sum(spatloc==all_vols(vols))==1; text = [volume ' sloc'];
                elseif sum(digloc==all_vols(vols))==1; text = [volume ' dig'];
                elseif sum(mdloc==all_vols(vols))==1; text = [volume ' md'];

                end
                %write next line to 'runs.list'
                fprintf(runslist,'%s\n',text);
            end
            fclose(runslist);
        end

        % Read in the runs.list for this session
        runnrs_tmp = [];
        fid = fopen([char(raw_path(session)), 'runs.list']);
        line = fgetl(fid);
        while ischar(line) && ~isempty(line)
            runnrs_tmp = [runnrs_tmp; line(1:3)];
            line = fgetl(fid);
        end
        fclose(fid);
        runnrs{session} = runnrs_tmp;

        % check if for this session the number of runs on my "runs.list" matches the number of volumes found in "invols"
        if size(invols{session},1) ~= size(runnrs{session},1)
            fprintf('\n  WARNING: Number of raw data volumes does not match the number of runs on your runs.list...\n  ...RETURNING...\n');
            return;
        end

    end

    if make_list_big
        if exist([subpath 'runs.list'],'file')
            unix(['rm ' subpath 'runs.list'])
        end
        % Read in all sessions' runs.list files again and save information into the one big text file
        fid_allrunslist = fopen([subpath, 'runs.list'],'w+'); % make a single text file where info on *all* runs (of all sessions) will be kept
        for session = all_sessions

            fin = fopen([char(raw_path(session)), 'runs.list'],'r');
            while 1
                inline = fgetl(fin);
                if ~ischar(inline)
                    break;
                else
                    new_line = [sprintf('%02d',session), '_', inline];
                    fprintf(fid_allrunslist,[new_line, '\n']);
                end
            end
            fclose(fin);

        end % session loop
        fclose(fid_allrunslist);
    end

    %% step 1: coregistration of all sessions to Doreti anatomy
    if WannaCoregister
        for session = preproc_sessions

            fprintf(['\n\n##### COREGISTRATION: SUBJECT ' num2str(sub) ', SESSION ' num2str(session) ' #####\n\n'])

            % Create a MCTemplate for this session       
            if exist([subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz'], 'file') % does a template already exist?
                MCTemplate = false;
                rsp = input('MCTemplate volume already exists. Overwrite? (y/n) ', 's');
                if strcmp(rsp, 'y')
                    MCTemplate = true;
                end
            else
                MCTemplate = true;
            end        
            if MCTemplate % if a MCtemplate should be made, create a new template or overwrite an existing template (in 1 of 2 possible ways)
                rsp = input('Do you want to use the 1st timepoint of the 1st functional run for your template? (y/n) ', 's'); % one way is to use the first volume of the first run in this session
                if strcmp(rsp, 'y')
                    disp(' ... making your Motion Correction Template... ');
                    [s,~] = unix(['fslroi ', char(raw_path(session)), invols{session}(1).name, ' ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz 0 1'], '-echo');
                    if s==1 % a non-zero value returned by unix means the command has failed
                        fprintf('\n\nERROR: MC template creation failed. returning.\n\n');
                        return;
                    end
                else
                    rsp = input('Do you want to use the mean of 1st functional run for your template? (y/n) ', 's'); % another way is to use the mean of all volumes in the first run in this session
                    if strcmp(rsp,'y')
                        [s,~] = unix(['fslmaths ', char(raw_path(session)), invols{session}(1).name, ' -Tmean ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz'], '-echo');
                        if s==1 % a non-zero value returned by unix means the command has failed
                            fprintf('\n\nERROR: MC template creation failed. returning.\n\n');
                            return;
                        end
                    else
                        fprintf('Currently no other options are supported by this function, returning.\n\n');
                        return;
                    end
                end
            end        

            % Manual registration between the MCTemplate from this session and
            % the anatomical volume from the retinotopy session
            ManReg = true;
            if exist([subpath 'registration_manual', sprintf('%02d',session), '.dat'], 'file') %does a manual registration already exists?
                rsp = input('A manual registration file already exists for this subject and session, do you want to overwrite? (y/n) ', 's');
                if strcmp(rsp,'n')
                    ManReg = false;
                end
            end
            if ManReg  %if I wanna do a manual registration, do a manual registration
                [s,w] = unix(['tkregister2 --s ' sub_doreti ' --mov ' subpath 'MCTemplate', sprintf('%02d',session), '.nii.gz --regheader-center --reg ' subpath 'registration_manual', sprintf('%02d',session), '.dat']);
                if s == 1 %a non-zero value returned by unix means the command has failed
                    fprintf(['\n\nERROR: Failed to start manual registration. \n\n' 'ERROR MESSAGE FROM UNIX WRAPPER:' w '\n\n']);
                    return
                end
            end

            % Automatic boundary based registration which will use the manual
            % registration as its input to do a superduper good finetuning 
            AutoReg = true;
            if exist([subpath 'registration_auto', sprintf('%02d',session), '.dat'], 'file') %does an automatic registration already exist?
                rsp = input('An automatic registration file already exists for this subject and session, do you want to overwrite? (y/n) ', 's');
                if strcmp(rsp,'n')
                    AutoReg = false;
                end
            end
            if AutoReg %if I wanna do an automatic registration, do an automatic registration
                fprintf('\n...Performing automatic registration using boundary based alignment, give it a minute...\n');
                unix(['bbregister --s ' sub_doreti ' --mov ' subpath 'MCTemplate', sprintf('%02d',session), '.nii.gz --reg ',...
                    subpath 'registration_auto', sprintf('%02d',session), '.dat --init-reg ' subpath 'registration_manual', ...
                    sprintf('%02d',session), '.dat --bold --o ' subpath 'bboutput', sprintf('%02d',session), '.nii.gz']);
                fprintf('Finished automatic registration!\n');
            end

            % Have the option to check the registration visually using tkregister
            % You will get an outline of the segmentation as well, this is because
            % FreeSurfer uses this segmentation to do the boundary based registration.
            check_auto = input('Do you want to visually inspect the result of the automatic registration? (y/n) ', 's');
            if strcmp(check_auto,'y')
                unix(['tkregisterfv --mov ' subpath 'MCTemplate', sprintf('%02d',session), '.nii.gz --reg ' subpath 'registration_auto', sprintf('%02d',session), '.dat --surfs']);
            end

            % Create the correct transformation matrices...
            % First --> create the correct "header transform" matrix
            if session == 1 % only need to do this once, because tkregister2 used the --regheader-center flag, thus all your header transforms will be identical across days
                CreateXFM = true;
                if exist([subpath 'regheadercenterMOD.mat'], 'file') % does a modified header file already exist?
                    rsp = input('A modified header registration file already exists for this subject, do you want to overwrite? (y/n) ', 's');
                    if strcmp(rsp,'n')
                        CreateXFM = false;
                    end
                end           
                if CreateXFM
                    % (1) create header-based transformation matrix files for session 1, based on a centered alignment between MCTemplate and orig.mgz
                    unix(['tkregister2 --s ', sub_doreti, ' --mov ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz --reg ', subpath, 'regheadercenter.dat --fslregout ', subpath, 'regheadercenter.mat --regheader-center --noedit']);
                    % (2) go in and check that this is indeed correct, i.e. that the brain is in the box (this depends on how the orig.mgz was acquired)
                    unix(['flirt -in ', subpath, 'bboutput', sprintf('%02d',session), '.nii.gz -ref ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz -applyxfm -init ', subpath, 'regheadercenter.mat -out ', subpath, 'TMP.nii.gz']);
                    rsp = 'c';
                    while ~strcmp(rsp,'h') %while not 'h'appy with how the brain is in the box
                        rsp = input(['\n --> Waiting for you to check ''TMP.nii.gz'' (in your subjects DataPreproc folder) with fslview... \n\n'...
                            ' Modify the ''regheadercenter.mat'' file if necessary and save as ''regheadercenterMOD'' \n'...
                            ' Remember: Change numbers (in steps of ~5) in the last column of the .mat to shift the brain \n'...
                            '           --> change number in top row to move brain sideways (probably never required) \n'...
                            '           --> add to middle row to move brain backward (subtract to move forward) \n'...
                            '           --> add to bottom row to move brain up (subtract to move down) \n\n'...
                            ' Press ''c'' to change, or ''h'' if you were happy with your new brain --> '],'s');
                        if strcmp(rsp,'c')
                            disp(' --> changing your box!')
                            unix(['flirt -in ', subpath, 'bboutput', sprintf('%02d',session), '.nii.gz -ref ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', subpath, 'TMP.nii.gz']);
                        end
                    end
                    unix(['rm ', subpath, 'TMP.nii.gz']); %clean up TMP volume
                end
            end

            % Second --> Transform my MCTemplate for this session (just to have it and to compare it in fsl view - the xfm template for session 1 == last TMP.nii.gz that was created).
            unix(['flirt -in ', subpath, 'bboutput', sprintf('%02d',session), '.nii.gz -ref ', subpath, 'MCTemplate', sprintf('%02d',session),...
                '.nii.gz -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', subpath, 'MCTemplateXFM', sprintf('%02d',session), '.nii.gz']);

            % Third --> I will need the output from FreeSurfer's bbregister (.dat) to be converted to FSL format (.mat)
            unix(['tkregister2 --s ', sub_doreti, ' --mov ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz --reg ', subpath,...
                'registration_auto', sprintf('%02d',session), '.dat --fslregout ', subpath, 'registration_auto', sprintf('%02d',session), '.mat --noedit']);

            % Finally, I will need the (MODIFIED) header info and the bbregister output matrix to be combined
            unix(['convert_xfm -omat ', subpath, 'concat', sprintf('%02d',session), '.mat -concat ', subpath, ...
                'regheadercenterMOD.mat ', subpath, 'registration_auto', sprintf('%02d',session), '.mat']);

        end %session loop (prior to other preproc step so that I can do all my manual work here before shit really hits the fan)
    end %coregistration loop


    %% step 2: transform all my functional data to common space
    if WannaTransform
        for session = preproc_sessions

            fprintf(['\n\n##### TRANSFORMATION: SUBJECT ' num2str(sub) ', SESSION ' num2str(session) ' #####\n\n'])

            % check that coregistration has actually finished for this session
            if ~exist([subpath, 'concat', sprintf('%02d',session), '.mat'], 'file')
                fprintf('\n  WARNING: It appears coregistration was incomplete...\n RETURNING\n\n');
                return;
            end

            % check if a transformation was already applied
            DoTrans = true;
            if exist([subpath, sprintf('%02d',session), '_REG_', runnrs{session}(end,:),'.nii.gz'],'file') %looks just for the last run
                rsp = input(['It appears nifti''s have already been transformed for subject ', num2str(sub), ', session ' num2str(session), '...',...
                    '\nDo you want to overwrite? (y/n) --> '],'s');
                if strcmp(rsp, 'n')
                    DoTrans = false;
                end
            end

            % loop over functional runs in this session, and transform the data
            % developers note: it is possible to skip all these transformations
            % by utilizing the -init flag in MCFLIRT to feed in the concat.mat 
            % as initial registration matrix during motion correction. Downside
            % is that this -init flag also applies the initial registration 
            % to the MCTemplate, which leads to some problems (ask RR if you're 
            % eager to learn more). So this is currently not recommended.
            if DoTrans
                for run = 1:size(runnrs{session},1)
                    fprintf(['\n\nApplying transformation to subject ' num2str(sub) ', session ' num2str(session) ', run ', num2str(run), ' (of ', num2str(size(runnrs{session},1)), ' runs)\n\n'])
                    unix(['flirt -in ', char(raw_path{session}), char(invols{session}(run).name), ' -ref ', subpath, 'MCTemplate', sprintf('%02d',session), '.nii.gz',...
                        ' -out ', subpath, sprintf('%02d',session), '_REG_', runnrs{session}(run,:), '.nii.gz -init ', subpath, 'concat', sprintf('%02d',session), '.mat -applyxfm']);
                end
            end

        end %session loop
    end %transformation loop


    %% step 3: motion correct all my functional data
    if WannaMC    

        % make folder to put first volume of each run in (enables visual check of MC in fslview)
        if ~exist([subpath 'FirstVols'],'dir')
            unix(['mkdir ', subpath, 'FirstVols', '/']);
        end   

        % do this by session, so you have the option to MC 1 session at a time
        for session = preproc_sessions 

            fprintf(['\n\n##### MOTION CORRECTION: SUBJECT ' num2str(sub) ', SESSION ' num2str(session) ' #####\n\n'])

            % check that transformations have actually finished for this session
            if ~exist([subpath, sprintf('%02d',session), '_REG_', runnrs{session}(end,:),'.nii.gz'],'file') %looks just for the last run
                fprintf('\n  WARNING: It appears transformations were incomplete...\n\n  RETURNING...\n\n');
                return;
            end

            % determine the volumes to go into motion correction
            clear invols % invols was previously defined as the raw data
            invols = dir([subpath, sprintf('%02d',session), '_REG_*nii.gz']); % here we will instead use the transformed data
            % MMH adding this 5/13/19 - check to make sure we're only looking
            % at volumes which are NOT motion corrected yet. Otherwise if we
            % have a crash in the middle of this script and try to restart it, 
            % we'd get an error.
            invols = invols(~contains({invols.name},'MC'));
            if size(runnrs{session},1)~=length(invols) % make sure I have the right number of invols because I'm double checking everything like stupid
                fprintf(['  WARNING: The number of runs in your runs.list does not match the number of transformed niti''s for session ', num2str(session), '...\n\n  RETURNING\n\n']);
                return;
            end


            % see if for this session moco was already performed or not
            DoMoco = true;
            if ~isempty(dir([subpath, sprintf('%02d',session), '_REG_MC*.nii.gz']))
                rsp = input(char(sprintf('Motion corrected volumes from session %d already exist, Overwrite? (y/n) --> ',session )),'s');
                if strcmp(rsp,'n')
                    DoMoco = false;
                end
            end

            % do the moco
            VolumeEnds = zeros(1,length(invols)+1);
            if DoMoco
                for func_idx = 1:length(invols) %for each of the to-be-preprocessed volumes in this session
                    fprintf(['\n\nProcessing moco for volume ', invols(func_idx).name, '...\n\n']);
                    % do moco
                    [s,~] = unix(['mcflirt -in ', subpath, invols(func_idx).name, ' -smooth 0 -stages 4 -dof 12 -reffile ', subpath,...
                        'MCTemplateXFM01.nii.gz -out ', subpath, sprintf('%02d',session), '_REG_MC_', runnrs{session}(func_idx,:), '.nii.gz -plots -mats'], '-echo');
                    if s==1
                        fprintf(['\n\nERROR: Motion Correction failed whilst processing volume ', invols(func_idx).name, '. returning.\n\n']);
                        return;
                    end
                    % read in the length of each par file for plotting purposes later
                    fid = fopen([subpath, sprintf('%02d',session), '_REG_MC_', runnrs{session}(func_idx,:), '.nii.gz.par'] , 'r'); %open the .par file where the translations and rotations are stored
                    tmp = fscanf(fid, '%f', [6 inf]); %how long (how many volumes) was this run
                    VolumeEnds(func_idx+1) = VolumeEnds(func_idx)+size(tmp,2);
                    fclose(fid);
                    % get first volume and put in "FirstVols" folder
                    if ~exist([subpath, 'FirstVols/', sprintf('%02d',session), '_', runnrs{session}(func_idx,:), '_Vol1.nii.gz'],'file')
                        unix(['fslroi ', subpath, sprintf('%02d',session), '_REG_MC_', runnrs{session}(func_idx,:), '.nii.gz ', subpath, 'FirstVols/', sprintf('%02d',session), '_', runnrs{session}(func_idx,:), '_Vol1.nii.gz.nii.gz 0 1'], '-echo');
                    end
                end
            end

            % evaluate MC outcomes for this session (save out plots)
            if exist([subpath, sprintf('%02d',session), 'EstimatedMotion.txt'], 'file'), delete([subpath, sprintf('%02d',session), 'EstimatedMotion.txt']); end
            %concat all translations & rotations and store in a separate text file
            unix(['cat ', subpath, sprintf('%02d',session), '*MC*.par >> ', subpath, sprintf('%02d',session), 'EstimatedMotion.txt']);
            %Read in estimated motion for this entire session across all runs (to use it for making mc plots)
            fid = fopen([subpath, sprintf('%02d',session), 'EstimatedMotion.txt']);
            EstMot = fscanf(fid, '%f', [6 inf])';
            fclose(fid);
            % --> roll pitch yaw
            h = figure('Position', [200 200 600 400], 'Visible', 'Off');
            plot(EstMot(:,1:3)/pi*360);
            v = axis; axis([0 VolumeEnds(end) v(3:4)]);
            set(gca, 'Xtick', VolumeEnds, 'Box', 'Off', 'Xgrid', 'On');
            xlabel('TR', 'FontSize', 16); ylabel('Rotation (deg.)', 'FontSize', 16);
            title('MCFLIRT Estimated Rotations', 'FontSize', 20, 'FontWeight', 'b');
            lh = legend({'Roll', 'Pitch', 'Yaw'}, 'FontSize', 16);
            set(lh, 'Box', 'Off'); set(gcf, 'PaperPositionMode', 'auto');
            print('-dtiff ', [subpath, sprintf('%02d',session), 'MC-RPY.tif']); close(h);
            % xyz
            h = figure('Position', [200 200 600 400], 'Visible', 'Off');
            plot(EstMot(:,4:6));
            v = axis; axis([0 VolumeEnds(end) v(3:4)]);
            set(gca, 'Xtick', VolumeEnds, 'Box', 'Off', 'Xgrid', 'On');
            xlabel('TR', 'FontSize', 16); ylabel('Translation (mm.)', 'FontSize', 16);
            title('MCFLIRT Estimated Translations', 'FontSize', 20, 'FontWeight', 'b');
            lh = legend({'X', 'Y', 'Z'}, 'FontSize', 16);
            set(lh, 'Box', 'Off'); set(gcf, 'PaperPositionMode', 'auto');
            print('-dtiff ', [subpath, sprintf('%02d',session), 'MC-XYZ.tif']); close(h);      

        end % session loop

        % combine all first volumes (of all runs and all sessions) into a single .nii.gz to check how well my MC worked out
        fprintf('\n Merging first volumes of all runs across all sessions for moco QC...\n');
        [s,~] = unix(['fslmerge -t ' subpath 'FirstVols/concat_volume.nii.gz ' subpath, 'FirstVols/*Vol1.nii.gz'],'-echo');
        if s==1 % a non-zero value returned by unix means the command has failed
            fprintf('\n\nERROR: Merging 1st volumes of all runs failed. returning.\n\n');
            return;
        end

    end % motion correction loop


    if WannaDetrend

        % do this by session, so you have the option to detrend 1 session at a time
        for session = preproc_sessions

            fprintf(['\n\n##### DETRENDING: SUBJECT ' num2str(sub) ', SESSION ' num2str(session) ' #####\n\n'])

            % check that motion correction has actually finished for this session
            if ~exist([subpath, sprintf('%02d',session), '_REG_MC_', runnrs{session}(end,:),'.nii.gz'],'file') % looks just for the last run
                fprintf('\n  WARNING: It appears motion correction was incomplete...\n\n  RETURNING...\n\n');
                return;
            end

            % determine the volumes to be detrended (assuming all other proprocessing has been done already)
            clear invols % invols was previously defined as the transformed data
            invols = dir([subpath, sprintf('%02d',session), '_REG_MC_*nii.gz']); % this time I take all the motion corrected volumes (of all sessions combined)
            HighPassCutoff = 50; % volumes

            % Do the detrending here
            DoDet = true;
            if length(dir([subpath, sprintf('%02d',session), '_REG_MC_DET*.nii.gz'])) == length(invols) % check if detrending was already performed and if so if you want to overwrite
                rsp = input('Detrending seems to have completed already. Overwrite? (y/n) ', 's');
                if strcmp(rsp, 'n')
                    DoDet = false;
                end
            end
            if DoDet %go ahead and detrend every volume collected (for this session)
                for func_idx = 1:length(invols)
                    fprintf(['\n\nProcessing volume ', invols(func_idx).name, '...\n\n']);
                    nmidx = invols(func_idx).name; %name of the output file
                    % Detrending (high pass filering) also removes the mean (the
                    % slowest drift in your data, if you will). If you want to add
                    % the mean back in you can use fslmaths to grab the mean first
                    % before you detrend (and add it back in during the detrend).
                    unix(['fslmaths ', subpath, invols(func_idx).name, ' -Tmean ', subpath, 'tmpMean']);
                    % actual detrending
                    [s,~] = unix(['fslmaths ', subpath, invols(func_idx).name, ' -bptf ', num2str(HighPassCutoff), ' -1 -add ', subpath, 'tmpMean ', subpath, nmidx(1:10), 'DET_', nmidx(11:end)], '-echo');
                    if s==1
                        fprintf(['\n\nERROR: Detrending failed whilst processing volume ', invols(func_idx).name, '. returning.\n\n']);
                        return;
                    end
                end

            end
        end % session
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% THE END! HAPPY CONTINUED ANALYSES TO YOU! %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end