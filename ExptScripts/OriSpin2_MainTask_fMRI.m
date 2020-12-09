% OriSpin2 Main Task 
% Spatial WM task - judge which side of the boundary disk the memory item 
% was on. 

% MMH 2/1/20
% updated by MMH 8/28/20 to fix a bug in which the preview boundary
% orientation was not being saved for random trials. 
Screen('Preference', 'SkipSyncTests', 1);    % added by MMH to make this code run on mac laptop

try
    %This should clear out the memory and stop sourcecode from being reprinted
    echo off
    clear 
    close all hidden
       
    debug = 0;  % debug mode, timing is sped up.
    p.TrainingSubject = 0;  % set to 0 to run mri experiment, 1 to train w laptop keys.

    % Get user info
    Subject = 'S98';
    p.Name = input('Initials? (default is temp) --> ','s'); if isempty(p.Name); p.Name = 'tmp'; end % collect subject initials
    SubNum = input('Subject number? (default is "0") --> ','s'); if isempty(SubNum); SubNum = '00'; end % collect subject number
    p.SubNumStr = sprintf('%02d', str2double(SubNum));
    p.SubNum = str2double(SubNum);
    if ~strcmp(Subject,['S' p.SubNumStr]); disp('Subject name doesn''t match name in script, please check and try again'); return; end

    SessNum = input('Session? (default is "1" -->','s');
    if isempty(SessNum)
        error('enter a valid session number')
    end
    p.SessNum = str2double(SessNum);
    
    p.RunsPerSess = 10;
    
    RunNum = input('Run? (default is "1") --> ','s');
    if isempty(RunNum)
        error('enter a valid run number')
    end
    p.RunNumGlobal = str2double(RunNum);    % note global is across this session only
   
    %% OPEN/INIT DATA FILE
    
    % save a copy of the currently running script here (in text format) so
    % we can look at it later if issues arise
    p.MyScript = fileread([mfilename('fullpath'),'.m']);

    datadir_local = [pwd filesep 'Data' filesep];
    p.TheDate = datestr(now,'yymmdd'); %Collect todays date (in p)
    p.TimeStamp = datestr(now,'HHMM'); %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
    
    if p.TrainingSubject==1
        p.fnsave = [datadir_local 'OriSpin2_TRAINING_S', p.SubNumStr, '_' p.TheDate '_MainTask.mat'];
    else
        p.fnsave = [datadir_local 'OriSpin2_S', p.SubNumStr, '_' p.TheDate '_MainTask.mat'];
    end
    
    if exist(p.fnsave,'file')
        load(p.fnsave);
    else            
        trialInfo = OriSpin2_CounterBal_MainTask(p.SubNum,0);
    end

    if exist('TheData','var')
        if p.RunNumGlobal~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.RunNumGlobal~=1
        error('No data exists yet for this subject, check your run number')
    end
       
    %% EXTRACT INFORMATION FROM THE TRIAL INFO STRUCTURE
    
    p.NumTrials = trialInfo.nTrialsBlock;
    p.rndseed_cbscript = trialInfo.rndseed;
    p.nTasks = trialInfo.nTasks;
    % extract the info from the trialInfo structure that we want to use now
    p.Task = squeeze(trialInfo.Task(1,p.RunNumGlobal,:));  
    p.TargPosition = squeeze(trialInfo.TargPosition(1,p.RunNumGlobal,:));  
    p.BoundPosition = squeeze(trialInfo.BoundPosition(1,p.RunNumGlobal,:));
    p.RandBoundPosition = squeeze(trialInfo.RandBoundPosition(1,p.RunNumGlobal,:));
    p.BonusTrial = squeeze(trialInfo.BonusTrial(1,p.RunNumGlobal,:));
    % converting the probe boundary to psychtoolbox orientation system
    p.ProbeOrientPTB = mod((-1)*(p.BoundPosition-90),360);
    p.RandProbeOrientPTB = mod((-1)*(p.RandBoundPosition-90), 360);
    % for the correct response - 1 always means left, 2 always right. 
    % depending on the session, left/right luminance mapping is swapped.
    p.CorrectResponse = squeeze(trialInfo.CorrectResponse(1,p.RunNumGlobal,:));
    
    p.ConditionStrs = {'Predictable','Random'};
    p.FingerStrs = {'Left Index','Right Index'};
      
    p.NumTrials = trialInfo.nTrialsBlock;
    
    % Decide which colors indicate which conditions for this subject
    if mod(p.SubNum,2)
        p.ConditionCueColors = {[14.5922, 79.0745, 150.6510],[245.0980, 36.8980, 17.6039]};
    else
        p.ConditionCueColors = {[245.0980, 36.8980, 17.6039],[14.5922, 79.0745, 150.6510]};

    end    
    
    % Luminance/finger mapping: odd sessions get swapped fingers.
    p.reverse_fingers = mod(p.SessNum,2);

    % this color is shown briefly before each target appears, to get ready
    p.PreTargCueColor = [0 121 0];
    %% TIMING INFORMATION 
    % take out the t structure here, has a lot of important fields
    t = trialInfo.t;
  
    if debug==1
        % make things super fast if you want to debug code
        t.BeginFixation = 1;
        t.EndFixation = 1;
        t.PostCueDelayLength = 4;
    end
    t.ITI = squeeze(trialInfo.ITI(1,p.RunNumGlobal,:))';
  
 
    %% SCREEN PARAMETERS
    Screens = Screen('Screens'); %look at available screens
    ScreenNumber = Screens(1); %pick first screen 
    p.ScreenSizePixels = Screen('Rect', ScreenNumber);
    tmprect = get(0, 'ScreenSize');
    computer_res = tmprect(3:4);
    if computer_res(1) ~= p.ScreenSizePixels(3) || computer_res(2) ~= p.ScreenSizePixels(4)
        Screen('CloseAll');clear screen;ShowCursor;
        disp('*** ATTENTION *** Yo screensizes ain''t matchin''')
    end
    CenterXPix = p.ScreenSizePixels(3)/2;
    CenterYPix = p.ScreenSizePixels(4)/2;
    ScreenHeight = 16; % in cm 
    ViewDistance = 49; % in cm
    % these measurements define all conversions from deg to pixels
    p.VisAngle = (2*atan2(ScreenHeight/2, ViewDistance))*(180/pi); % visual angle of the whole screen
    p.ppd = p.ScreenSizePixels(4)/p.VisAngle; % pixels per degree visual angle
    p.MyGrey = 128;
    gammacorrect = false;

    
    %% KEYS
    KbName('UnifyKeyNames')

    % keys for scanner button box
    % Order from left to right is [b,y,g,r]
    p.keys=[KbName('b'),KbName('r')];
    if p.TrainingSubject==1
        p.keys = [KbName('g'),KbName('h')];
    end
    p.escape = KbName('escape');
    p.space = KbName('space');
    p.start = KbName('t');

    %% WINDOW SETUP
    AssertOpenGL;
    PsychJavaTrouble;
    [w] = Screen('OpenWindow',ScreenNumber, p.MyGrey, [], [],2);
    [p.ifi] = Screen('GetFlipInterval',w);
    HideCursor;
    OriginalCLUT = Screen('LoadClut', w);   % need this later on
    white=WhiteIndex(ScreenNumber);
    black=BlackIndex(ScreenNumber);
    

    %% STIMULUS/DISPLAY PARAMETERS
    % information about the target (memory) dot
    p.TargDistancePix = 7*p.ppd; % distance from fixation in pixels
    p.TargColor = white;
    p.TargRadiusPix = 0.15*p.ppd; % size of dot in pixels
    
    % information about fixation pt and aperture
    p.ApertureRadiusPix = 0.4*p.ppd; 
    p.FixRadiusPix = 0.2*p.ppd; % size of dot in pixels
    p.FixColor = black;
    %[left, top, right, bottom]
    p.FixRect = [CenterXPix-p.FixRadiusPix,...
        CenterYPix-p.FixRadiusPix,...
        CenterXPix+p.FixRadiusPix,...
        CenterYPix+p.FixRadiusPix];
    p.ApertureRect = [CenterXPix-p.ApertureRadiusPix,...
        CenterYPix-p.ApertureRadiusPix,...
        CenterXPix+p.ApertureRadiusPix,...
        CenterYPix+p.ApertureRadiusPix];

    % information about the instruction screen and cue dots
    space = 1;  % spacing between lines of text in degrees
    %[left, top, right, bottom]
    p.CueStartRects = repmat(p.FixRect, p.nTasks, 1) + ...
        [-space*p.ppd, -2*space*p.ppd, -space*p.ppd, -2*space*p.ppd;...
        -space*p.ppd, -space*p.ppd, -space*p.ppd, -space*p.ppd];
    %[x,y] in pixels
    p.CueStartTextPos = [p.CueStartRects(:,3)+space*p.ppd, p.CueStartRects(:,2)+0.35*p.ppd];
    
    p.ProbeStartRects = repmat(p.FixRect, 2, 1) + ...
        [-space*p.ppd, 1.5*space*p.ppd, -space*p.ppd, 1.5*space*p.ppd;...
        -space*p.ppd, 2.5*space*p.ppd, -space*p.ppd, 2.5*space*p.ppd];
        
    p.ProbeStartTextPos = [p.ProbeStartRects(:,3)+space*p.ppd, p.ProbeStartRects(:,2)+0.35*p.ppd];
        
    % information about the probe disk 
    p.ProbeRadiusDeg = 10;   
    p.ProbeRadiusPix = round(p.ProbeRadiusDeg*p.ppd); % entire size of the probe patch
    p.Smooth_size = round(0.5*p.ppd); %size of fspecial smoothing kernel
    p.Smooth_sd = round(0.5*p.ppd); %smoothing kernel sd    
    p.ProbeRadiusSmall = p.ProbeRadiusPix-p.Smooth_size/2; % size of the probe stimulus before it gets smoothed
    p.ApertureRadiusBig  = p.ApertureRadiusPix+p.Smooth_size/2;  % size of the aperture stimulus before it gets smoothed
    p.BoundContrast = 0.1;
    p.BoundSideColors = {(p.MyGrey-p.BoundContrast*255)*[1,1,1], ...
        (p.MyGrey+p.BoundContrast*255)*[1,1,1]};

    if p.reverse_fingers
        % want each color to mean opposite finger
        p.BoundSideColors = p.BoundSideColors([2,1]);
    end
    %% MAKE STIMULI 

    % define exact positions to draw memory dots at
    % making y negative here because psychtoolbox draws up as down
    p.TargPosX = p.TargDistancePix*cosd(p.TargPosition);
    p.TargPosY = -p.TargDistancePix*sind(p.TargPosition);
    %[-x,-y,x,y]
    p.TargRectsPix = repmat([CenterXPix, CenterYPix, CenterXPix, CenterYPix],p.NumTrials,1) +...
        repmat([-p.TargRadiusPix,-p.TargRadiusPix, p.TargRadiusPix, p.TargRadiusPix],p.NumTrials,1) +...
        [p.TargPosX,p.TargPosY,p.TargPosX,p.TargPosY]; 
    
    % define rectangle for where to show the disk (boundary)
    ProbeDiskRect = [CenterXPix-p.ProbeRadiusPix CenterYPix-p.ProbeRadiusPix CenterXPix+p.ProbeRadiusPix CenterYPix+p.ProbeRadiusPix];
    % start with a meshgrid
    X=-p.ProbeRadiusPix+.5:1:p.ProbeRadiusPix-.5; Y=-p.ProbeRadiusPix+.5:1:p.ProbeRadiusPix-.5;
    [x,y] = meshgrid(X,Y);
    r = sqrt(x.^2 + y.^2);

    % create a donut stimulus cut into two halves
    smoothed_donut = zeros(size(x));
    smoothed_donut(r<=p.ProbeRadiusSmall & r>=p.ApertureRadiusBig) = 1;
    smoothed_donut = filter2(fspecial('gaussian', p.Smooth_size, p.Smooth_sd), smoothed_donut);
    
    halves_stim = zeros(size(x));
    halves_stim(x<0) = p.BoundSideColors{1}(1);
    halves_stim(x>0) = p.BoundSideColors{2}(1);
    
    probe_stim = halves_stim.*smoothed_donut;
    probe_stim = probe_stim+p.MyGrey*(1-smoothed_donut);

    p.probe_image = repmat(probe_stim,1,1,3);
    
    ProbeTexture = Screen('MakeTexture', w, p.probe_image);
    
    %% Allocate arrays to store trial info
    data.Response = nan(p.NumTrials,1); % this is the button they pressed
    
    t.stim_flips = nan(p.NumTrials,2); % stim on/stim off
    t.bound_flips = nan(p.NumTrials,2);  % bound on/bound off
    t.respTimeFromStart = nan(p.NumTrials,1);
       
     %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    Screen(w,'TextFont','Helvetica');

    % draw reminders of which cue color is which
    for ii=1:p.nTasks
        Screen('FillOval', w, p.ConditionCueColors{ii}, p.CueStartRects(ii,:));
        DrawFormattedText(w, p.ConditionStrs{ii}, p.CueStartTextPos(ii,1), p.CueStartTextPos(ii,2), white);    
    end
    % draw reminders of which finger is which
    for ii=1:2
        Screen('FillOval', w, p.BoundSideColors{ii}, p.ProbeStartRects(ii,:));
        DrawFormattedText(w, p.FingerStrs{ii}, p.ProbeStartTextPos(ii,1), p.ProbeStartTextPos(ii,2), white);    
    end

    Screen('Flip', w);
%     ListenChar(2) % added by MMH to make this run on mac laptop
    resp=0;
    % wait for a space bar press to start
    while resp==0
        [resp, timeStamp] = checkForResp([p.start,p.space],p.escape);

        if resp==-1; escaperesponse(OriginalCLUT); end    
    end
    t.StartTime = GetSecs;
    
    KbReleaseWait();
    
    %% Fixation period before starting the stimuli (for scanner, this is the 12.8 seconds thing)
    FlushEvents('keyDown');
    Screen('Flip', w);
    ListenChar(2)
    
%     t.StartTime = GetSecs; %Get the starttime of the experiment in seconds
    GlobalTimer = 0; %this timer keeps track of all the timing in the experiment. TOTAL timing.
    TimeUpdate = t.StartTime; %what time is it now?
    % presentt begin fixation
    Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
    Screen('FillOval', w, p.FixColor, p.FixRect);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    %TIMING!:
    GlobalTimer = GlobalTimer + t.BeginFixation;
    TimePassed = (GetSecs-TimeUpdate); %Flush the time the previous event took
    while (TimePassed<t.BeginFixation) %For as long as the cues are on the screen...
        TimePassed = (GetSecs-TimeUpdate);%And determine exactly how much time has passed since the start of the expt.       
        [resp, ~] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
    end
    TimeUpdate = TimeUpdate + t.BeginFixation;

    %% start trial loop

    for tt = 1:p.NumTrials

        t.TrialStartTime(tt) = GlobalTimer; %Get the starttime of each single trial (relative to experiment start)
        
        %% Pre-target fixation cue
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);       
        Screen('FillOval', w, p.PreTargCueColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);

        %TIMING!:
        GlobalTimer = GlobalTimer + t.PreTargetCueTime;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.PreTargetCueTime)
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.PreTargetCueTime; %Update Matlab on what time it is.
   
        
        %% Target (spatial memory item)

        Screen('FillOval', w, p.TargColor, p.TargRectsPix(tt,:))
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);       
        Screen('FillOval', w, p.PreTargCueColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        t.stim_flips(tt,1) = GetSecs;

        %TIMING!:
        GlobalTimer = GlobalTimer + t.TargetTime;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.TargetTime)
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.TargetTime; %Update Matlab on what time it is.

      
        %% delay 1 (short) 
        
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        t.stim_flips(tt,2) = GetSecs;
        %TIMING!:
        GlobalTimer = GlobalTimer + t.PreCueDelayLength;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.PreCueDelayLength) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.PreCueDelayLength; %Update Matlab on what time it is.
       
        %% Condition Cue (Predictable/Random)
        
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.ConditionCueColors{p.Task(tt)}, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
      
        %TIMING!:
        GlobalTimer = GlobalTimer + t.CueLength;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.CueLength) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.CueLength; %Update Matlab on what time it is.
       
        %% show a boundary disk - random or informative
        
        if p.Task(tt)==1
            % show the real orient
            orient_to_draw = p.ProbeOrientPTB(tt);
        else
             % show a random orient
            orient_to_draw = p.RandProbeOrientPTB(tt);
%             orient_to_draw = datasample(trialInfo.binListBound(:),1);
        end
        
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillRect', w, p.MyGrey);
        Screen('DrawTexture', w, ProbeTexture, [], ProbeDiskRect, orient_to_draw, 0);
        Screen('FillOval', w, p.MyGrey, p.ApertureRect);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        
        GlobalTimer = GlobalTimer + t.BoundPreviewLength;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.BoundPreviewLength) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.BoundPreviewLength; %Update Matlab on what time it is.
       
        %% delay 2 (long)
        
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);

        %TIMING!:
        GlobalTimer = GlobalTimer + t.PostCueDelayLength;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.PostCueDelayLength) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
        end
        TimeUpdate = TimeUpdate + t.PostCueDelayLength; %Update Matlab on what time it is.
       
         
        %% boundary disk - response period now

        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillRect', w, p.MyGrey);
        Screen('DrawTexture', w, ProbeTexture, [], ProbeDiskRect, p.ProbeOrientPTB(tt), 0);
        Screen('FillOval', w, p.MyGrey, p.ApertureRect);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        % first flip of the disk is the beginning of their response period
        t.bound_flips(tt,1) = GetSecs;
        %TIMING!:
        GlobalTimer = GlobalTimer + t.BoundFinalLength;
        TimePassed = (GetSecs-TimeUpdate); 
        respYet = 0;
        while (TimePassed<t.BoundFinalLength)
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end
            if ~respYet
                if ismember(resp,p.keys)
                    % response! 
                    respYet=1;
                    t.respTimeFromStart(tt) = GetSecs - t.bound_flips(tt,1);
                    data.Response(tt) = find(p.keys==resp);           
                end
            end
        end
        TimeUpdate = TimeUpdate + t.BoundFinalLength; %Update Matlab on what time it is.
        
        %% ITI  

        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);

        % record the dial offset time here
        t.bound_flips(tt,2) = GetSecs;
        
        %TIMING!:
        GlobalTimer = GlobalTimer + t.ITI(tt);
        TimePassed = GetSecs-TimeUpdate; %Flush time passed.
        while (TimePassed<t.ITI(tt)) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end
            % continue checking response into the ITI if they haven't
            % given a response yet, and until t.MaxRespTime is reached.
            if ~respYet && ((GetSecs - t.bound_flips(tt,1))<t.MaxRespTime)
                if ismember(resp,p.keys)
                    % response! 
                    respYet=1;
                    t.respTimeFromStart(tt) = GetSecs - t.bound_flips(tt,1);
                    data.Response(tt) = find(p.keys==resp);           
                end
            end
        end
        TimeUpdate = TimeUpdate + t.ITI(tt); %Update Matlab on what time it is.
     
        
    end 
    %% finish experiment 
       
    % final fixation:
    Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
    Screen('FillOval', w, p.FixColor, p.FixRect);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    GlobalTimer = GlobalTimer + t.EndFixation;
    TimePassed = GetSecs-TimeUpdate;
    while (TimePassed<t.EndFixation) 
         TimePassed = GetSecs-TimeUpdate; 
        [resp, ~] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
    end
    
    t.EndTime = GetSecs; %Get endtime of the experiment in seconds
    t.TotalExpTime = (t.EndTime-t.StartTime); %Gets the duration of the total run.
    t.TotalExpTimeMins = t.TotalExpTime/60; %TOTAL exp time in mins including begin and end fixation.
    t.GlobalTimer = GlobalTimer; %Spits out the exp time in secs excluding begin and end fixation.

     %----------------------------------------------------------------------
    %LOOK AT BEHAVIORAL PERFORMANCE---------------------------------------
    %----------------------------------------------------------------------
    
    % in case it was aborted early, will just look at perf so far
    inds2use = ~isnan(t.bound_flips(:,2));
    data.Accuracy = mean(data.Response(inds2use)==p.CorrectResponse(inds2use));
    
    % also get bonus point accuracy - how did they do on the most hard
    % trials?
    inds_bonus = ~isnan(t.bound_flips(:,2)) & p.BonusTrial==1;
    data.NumBonus = sum(data.Response(inds_bonus)==p.CorrectResponse(inds_bonus));
    
    fprintf('\nCompleted block %d!\n',p.RunNumGlobal);
    fprintf('Accuracy: %.2f percent\n',data.Accuracy*100);
    fprintf('Bonus points earned: %d\n',data.NumBonus);
    
    InstrText = ['Block finished!' '\n\n'...
                'Average accuracy: ' sprintf('%.2f',data.Accuracy*100) ' percent \n\n'...
                'Bonus points earned: ' sprintf('%d',data.NumBonus) '\n\n'...
                'Press space to exit'];

    DrawFormattedText(w, InstrText, 'center', 'center', white);
    % put up a message to wait
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
     
    %----------------------------------------------------------------------
    %SAVE OUT THE DATA-----------------------------------------------------
    %----------------------------------------------------------------------
%     if exist(p.fnsave,'file')
%        load(p.fnsave);
%     end
    %First I make a list of variables to save:
    TheData(p.RunNumGlobal).t = t;
    TheData(p.RunNumGlobal).p = p;
    TheData(p.RunNumGlobal).data = data;
    
    save(p.fnsave,'TheData','trialInfo');
    
    resp=0;
    % wait for a space bar press to exit
    while resp~=p.space
        [resp, timeStamp] = checkForResp(p.space, p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
    end

    KbReleaseWait();
    
    
    %----------------------------------------------------------------------
    %WINDOW CLEANUP--------------------------------------------------------
    %----------------------------------------------------------------------
    %This closes all visible and invisible screens and puts the mouse cursor
    %back on the screen
    Screen('CloseAll');
    if exist('OriginalCLUT','var')
        if exist('ScreenNumber','var')
            Screen('LoadCLUT', ScreenNumber, OriginalCLUT);
        else
            Screen('LoadCLUT', 0, OriginalCLUT);
        end
    end
    clear screen
    ListenChar(1);
    ShowCursor;
    
%----------------------------------------------------------------------
%TRY CATCH STUFF-------------------------------------------------------
%----------------------------------------------------------------------
catch err%If an error occurred in the "try" block, this code is executed
   
    if exist('OriginalCLUT','var') && ~isempty(OriginalCLUT)
        if exist('ScreenNumber','var')
            Screen('LoadCLUT', ScreenNumber, OriginalCLUT);
        else
            Screen('LoadCLUT', 0, OriginalCLUT);
        end
    end
    Screen('CloseAll');                
    ShowCursor;
    if IsWin
        ShowHideWinTaskbarMex;     
    end
    ListenChar(1)
    rethrow(err)

end
