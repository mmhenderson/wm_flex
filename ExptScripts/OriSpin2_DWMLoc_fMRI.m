% OriSpin2 Digit WM Localizer Task
% Delayed button pressing task - see a color which indicates which finger
% to press, withhold the press until after a delay period. 

% MMH 8/11/20

Screen('Preference', 'SkipSyncTests', 1);% added by MMH to make this code run on mac laptop

try
    %This should clear out the memory and stop sourcecode from being reprinted
    echo off
    clear 
    close all hidden
       
    debug = 0; % debug mode, timing is sped up.
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
        p.fnsave = [datadir_local 'OriSpin2_TRAINING_S', p.SubNumStr, '_' p.TheDate '_DWMLoc.mat'];
    else
        p.fnsave = [datadir_local 'OriSpin2_S', p.SubNumStr, '_' p.TheDate '_DWMLoc.mat'];
    end
    
    if exist(p.fnsave,'file')
        load(p.fnsave);        
    end
    
    % run this function on every single run, does all the 
    % randomization/counter-balancing. 
    % The trialInfo struct for each run gets saved in TheData(p.RunNumGlobal).trialInfo            
    trialInfo = OriSpin2_CounterBal_DWMLoc(p.SubNum);
     
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

    % extract the info from the trialInfo structure that we want to use now
    p.WhichFinger = squeeze(trialInfo.WhichFinger(1,:))'; 
    if mod(p.SubNum,2)  
        % odd subject, color mappings go [1,0]
        p.ColorOrder = mod(p.RunNumGlobal,2);
    else        
        % even subject, color mappings go [0,1]
        p.ColorOrder = 1-mod(p.RunNumGlobal,2);
    end  
    p.NumTrials = trialInfo.nTrialsBlock;

    % decide what the color to finger mapping is for this run
    if p.ColorOrder==1
        p.LeftRightColor = [1.5922 79.0745 155; 200 0 226.6; 80 80 80]; % left = blue, right = magenta
    else 
        p.LeftRightColor = [200 0 226.6; 1.5922 79.0745 155; 80 80 80]; % left = magenta, right = blue
    end
    %% TIMING INFORMATION 
    % take out the t structure here, has a lot of important fields
    t = trialInfo.t;
  
    if debug==1
        % make things super fast if you want to debug code
        t.BeginFixation = 1;
        t.EndFixation = 1;
        t.DelayLength = 4;
    end
    t.ITI = squeeze(trialInfo.ITI(1,:))';
  
 
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

    % what color the dot becomes when it is time to respond
    p.ResponseColor = white;
    
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

    %% Allocate arrays to store trial info
    data.Response = nan(p.NumTrials,1); % this is button they pressed
    data.RespTooEarly = zeros(p.NumTrials,1); 
    % this is 1 if they press a button any time after color change and before go cue.
    t.respTimeFromStart = nan(p.NumTrials,1);
    
    t.stim_flips = nan(p.NumTrials,2); % stim on/stim off
    t.resp_flips = nan(p.NumTrials,2); % response probe on/off
     
     %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    Screen(w,'TextFont','Helvetica');

    % draw instructions text - show finger/color mapping here.
    space=1;
    %[left, top, right, bottom]
    instr_dot_rects = repmat(p.FixRect, 2, 1) + ...
        [-space*p.ppd, -2*space*p.ppd, -space*p.ppd, -2*space*p.ppd;...
        -space*p.ppd, -space*p.ppd, -space*p.ppd, -space*p.ppd];
    %[x,y] in pixels
    instr_text_pos = [instr_dot_rects(:,3)+space*p.ppd, mean(instr_dot_rects(:,[2]),2)];
    
    Screen('FillOval', w, p.LeftRightColor(1,:), instr_dot_rects(1,:));
    Screen('FillOval', w, p.LeftRightColor(2,:), instr_dot_rects(2,:));
    
    Screen('DrawText',w, 'Left', instr_text_pos(1,1), instr_text_pos(1,2));
    Screen('DrawText',w, 'Right', instr_text_pos(2,1), instr_text_pos(2,2));
     
    Screen('FillOval', w, p.FixColor, p.FixRect);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
%     ListenChar(2);
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
        
        
        %% Color cue

        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);       
        Screen('FillOval', w, p.LeftRightColor(p.WhichFinger(tt),:), p.FixRect);
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
            if ismember(resp,p.keys)
                 % this is a too-early response!  mark it
                 data.RespTooEarly(tt) = 1;
            end
        end
        TimeUpdate = TimeUpdate + t.TargetTime; %Update Matlab on what time it is.

      
        %% delay period 
        
        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        t.stim_flips(tt,2) = GetSecs;
        %TIMING!:
        GlobalTimer = GlobalTimer + t.DelayLength;
        TimePassed = (GetSecs-TimeUpdate); 
        while (TimePassed<t.DelayLength) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end 
            if ismember(resp,p.keys)
                 % this is a too-early response!  mark it
                 data.RespTooEarly(tt) = 1;
            end
        end
        TimeUpdate = TimeUpdate + t.DelayLength; %Update Matlab on what time it is.
       
         
        %% response period - press L or R button now

        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.MyGrey, p.ApertureRect);
        Screen('FillOval', w, p.ResponseColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        % first flip of the disk is the beginning of their response period
        t.resp_flips(tt,1) = GetSecs;
        %TIMING!:
        GlobalTimer = GlobalTimer + t.MaxRespTime;
        TimePassed = (GetSecs-TimeUpdate); 
        respYet = 0;
        while (TimePassed<t.MaxRespTime)
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end
            if ~respYet
                if ismember(resp,p.keys)
                    % response! 
                    respYet=1;
                    t.respTimeFromStart(tt) = GetSecs - t.resp_flips(tt,1);
                    data.Response(tt) = find(p.keys==resp);           
                end
            end
        end
        TimeUpdate = TimeUpdate + t.MaxRespTime; %Update Matlab on what time it is.
        
        %% ITI  

        Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
        Screen('FillOval', w, p.FixColor, p.FixRect);
        Screen('DrawingFinished', w);
        Screen('Flip', w);

        % record the response period offset time here
        t.resp_flips(tt,2) = GetSecs;
        
        %TIMING!:
        GlobalTimer = GlobalTimer + t.ITI(tt);
        TimePassed = GetSecs-TimeUpdate; %Flush time passed.
        while (TimePassed<t.ITI(tt)) 
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, ~] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end           
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
    
    % look at performance on all trials that were shown
    inds2use = ~isnan(t.resp_flips(:,2));
    correct = data.Response(inds2use)==p.WhichFinger(inds2use);
    % if they responded early, mark the trial as incorrect
    correct(data.RespTooEarly(inds2use)==1) = 0;
    data.Accuracy = mean(correct);
    
    fprintf('\nCompleted block %d!\n',p.RunNumGlobal);
    fprintf('Accuracy: %.2f percent\n',data.Accuracy*100);
    
    InstrText = ['Block finished!' '\n\n'...
                'Accuracy: ' sprintf('%.2f percent',data.Accuracy*100) '\n\n'...               
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
    TheData(p.RunNumGlobal).trialInfo = trialInfo;
    
    save(p.fnsave,'TheData');
    
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


