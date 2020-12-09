%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!! This is a DIGIT LOCALIZER !!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% Random probes to the left and the right index finger are presented in a
% fast-event related design (3/5/8s iti). Colors (blue/magenta) indicate
% the finger to be used. Subject is instructed to press as soon as they see 
% the color change. The subject must maintain fixation throughout.
%
% Response mapping is counterbalanced between pairs of runs, make sure to
% do an even number within each session.
%
% CURRENTLY RUNS FOR 5:19 (319 seconds)
% TR is 800 ms (0.8 s)
% Must use 395 TR's for this run
%
% originally written by RR, Jul 2018.
% modified by MMH, Jan 2019.

%----------------------------------------------------------------------
%PREPARE AND COLLECT BASIC INFO----------------------------------------
%----------------------------------------------------------------------
echo off
clear
close all hidden

% Get path info
expdir = pwd;
datadir = 'Data';
GeneralUseScripts = '/home/serencesadmin/Documents/MATLAB/Rosanne/GeneralUseScripts/';
addpath(GeneralUseScripts); % add my general use scripts to the path.

% Get user info
Subject = 'S99';
p.Name = input('Initials? (default is temp) --> ','s'); if isempty(p.Name); p.Name = 'tmp'; end % collect subject initials
SubNum = input('Subject number? (default is "0") --> ','s'); if isempty(SubNum); SubNum = '00'; end % collect subject number
p.SubNum = sprintf('%02d', str2double(SubNum));
if ~strcmp(Subject,['S' p.SubNum]); disp('Subject name doesn''t match name in script, please check and try again'); return; end
p.MyScript = fileread([mfilename('fullpath') '.m']);

% Get timing info
t.MySeed = sum(100*clock); % seeds the random number generator based on the current time
rng(t.MySeed,'twister'); % sets random seed
t.TheDate = datestr(now,'yymmdd'); % collect todays date (in t.)
t.TimeStamp = datestr(now,'HHMM'); % timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)



%----------------------------------------------------------------------
%SCREEN PARAMETERS-----------------------------------------------------
%----------------------------------------------------------------------
Screens = Screen('Screens'); % look at available screens
ScreenNr = Screens(1); % pick first screen
p.ScreenSizePixels = Screen('Rect', ScreenNr); % size of the screen I will draw to
tmprect = get(0, 'ScreenSize'); computer_res = tmprect(3:4); % size of my first screen
if computer_res(1) ~= p.ScreenSizePixels(3) || computer_res(2) ~= p.ScreenSizePixels(4)
    Screen('CloseAll'); clear screen; ShowCursor;
    disp('*** ATTENTION! *** Yo screensizes ain''t matchin''')
end
CenterX = p.ScreenSizePixels(3)/2;
CenterY = p.ScreenSizePixels(4)/2;
ScreenHeight = 16; % in cm 
ViewDistance = 49; % in cm
%TODO: check these measurements!!! These numbers set the conversion
%from pixels to degrees, so if we want the degrees to stay constant
%across subs we must measure every time.
p.VisAngle = (2*atan2(ScreenHeight/2, ViewDistance))*(180/pi); % visual angle of the whole screen
p.ppd = p.ScreenSizePixels(4)/p.VisAngle; % pixels per degree of visual angle

p.MyGrey = 128;
black = BlackIndex(ScreenNr); white = WhiteIndex(ScreenNr);
gammacorrect = false;

% scanner buttons are : [b,y,g,r] from left to right
% we only care about the furthest left and furthest right.
p.keys=[KbName('b'),KbName('y'),KbName('g'),KbName('r')];

p.escape = KbName('escape');
p.space = KbName('space');
p.start = KbName('t');

%----------------------------------------------------------------------
%OPEN/INIT DATA FILES--------------------------------------------------
%----------------------------------------------------------------------
cd(datadir); % and go there to fetch stuff
if exist(['OriSpin_S', p.SubNum, '_DigitLocalizer.mat'],'file')
    load(['OriSpin_S', p.SubNum, '_DigitLocalizer.mat']);
    runnumber = length(TheData)+1; % set the number of the current run
    p.TrialNumGlobal = TheData(end).p.TrialNumGlobal;
    % switch to the opposite color mapping from the previous run.
    p.LeftRightColor = [flipud(TheData(end).p.LeftRightColor(1:2,:)); 80 80 80];
   
else
    runnumber = 1; % if no data file exists this must be the first run
    p.TrialNumGlobal = 0;
    my_flip = CoinFlip(1,.5);
    % randomly select which color mapping to start with.
    if rem(my_flip,2) % if uneven coinflip
        p.LeftRightColor = [1.5922 79.0745 155; 200 0 226.6; 80 80 80]; % left = blue, right = magenta
    else % if even coinflip
        p.LeftRightColor = [200 0 226.6; 1.5922 79.0745 155; 80 80 80]; % left = magenta, right = blue
    end
end
cd(expdir)




%----------------------------------------------------------------------
%DEFINE MAIN PARAMETERS------------------------------------------------
%----------------------------------------------------------------------

p.NumTrials = 60; % must be even number

% Stimulus params
p.OuterFixationRadius = .2*p.ppd; % outer dot radius
p.InnerFixationRadius = 0*p.ppd; % inner dot radius
p.FixColor = black;

% Timing params
t.BeginFixation =  12.8 + 1; % 16 TRs need to be extra for 'calibration' scans (16trs * .8ms);
t.InstructionTime = 4; % instructions are shown in this earliest portion of the BeginFixation
t.CueTime = 1; % duration of the finger-press cue
t.ITIrange = [2,6];
t.ITIs = linspace(t.ITIrange(1), t.ITIrange(2), p.NumTrials);

% randomly assign which finger will be cued on each trial
% 1=Left, 2=Right
p.WhichFinger = repmat([1;2],p.NumTrials/2,1);
p.WhichFinger = p.WhichFinger(randperm(p.NumTrials));

t.EndFixation = 5;

t.MeantToBeTime = t.BeginFixation + t.CueTime*p.NumTrials + sum(t.ITIs)+ t.EndFixation;

%----------------------------------------------------------------------
%WINDOW SETUP & GAMMA CORRECTION---------------------------------------
%----------------------------------------------------------------------
AssertOpenGL;
PsychJavaTrouble;
[w] = Screen('OpenWindow',ScreenNr, p.MyGrey,[],[],2);
t.ifi = Screen('GetFlipInterval',w);
OriginalCLUT = Screen('LoadClut', w);
if gammacorrect
    OriginalCLUT = Screen('LoadClut', w);
    MyCLUT = zeros(256,3); MinLum = 0; MaxLum = 1;
    CalibrationFile = 'calib_07-Nov-2016.mat';
    [gamInverse,dacsize] = LoadCalibrationFileRR(CalibrationFile, expdir, GeneralUseScripts);
    LumSteps = linspace(MinLum, MaxLum, 256)';
    MyCLUT(:,:) = repmat(LumSteps, [1 3]);
    MyCLUT = round(map2map(MyCLUT, repmat(gamInverse(:,4),[1 3]))); %Now the screen output luminance per pixel is linear!
    Screen('LoadCLUT', w, MyCLUT);
    clear CalibrationFile gamInverse
end
HideCursor;



%----------------------------------------------------------------------
%MISC AND PREALLOCATE STUFF--------------------------------------------
%----------------------------------------------------------------------
data.RespOnTime = zeros(p.NumTrials,1); % true if they respond anytime during the 1 second interval
data.Response = nan(p.NumTrials,1); % a 1 or a 4 (hopefully)
data.RespTimeRaw = nan(p.NumTrials,1); 
data.RespTimeFromOnset = nan(p.NumTrials,1); % reaction time
data.Correct = zeros(p.NumTrials,1); % correct button, correct time.

% record if they do anything wrong
data.RespLate = zeros(p.NumTrials,1);
data.WrongButton = zeros(p.NumTrials,1);    

% record onset/offset time of the cues
t.FlipTimes = nan(p.NumTrials,2);

%----------------------------------------------------------------------
%WELCOME MESSAGE & WAIT FOR TRIGGER------------------------------------
%----------------------------------------------------------------------
Screen(w,'TextSize',20);
Screen('DrawText',w, 'waiting for scanner to initiate', p.ScreenSizePixels(1)+10, p.ScreenSizePixels(2)+10, white);

%Present beginning of the begin fixation (with instructions)
Screen('FillOval', w, p.LeftRightColor(1,:), [CenterX-p.OuterFixationRadius-40 CenterY-p.OuterFixationRadius+25 CenterX+p.OuterFixationRadius-40 CenterY+p.OuterFixationRadius+25])
Screen('FillOval', w, p.LeftRightColor(2,:), [CenterX-p.OuterFixationRadius+20 CenterY-p.OuterFixationRadius+25 CenterX+p.OuterFixationRadius+20 CenterY+p.OuterFixationRadius+25])
Screen('DrawText',w, 'left', CenterX+p.OuterFixationRadius-38, CenterY-p.OuterFixationRadius*1.5+27, 20);
Screen('DrawText',w, 'right', CenterX+p.OuterFixationRadius+25, CenterY-p.OuterFixationRadius*1.5+27, 20);

Screen('Flip', w);
FlushEvents('keyDown'); %First discard all characters from the Event Manager queue.
ListenChar(2);
% Sittin' here, waitin' on my trigger...
resp=0;
while resp==0
    [resp, timestamp] = checkForResp([p.start,p.space],p.escape);
    if resp==-1; escaperesponse(OriginalCLUT); end 
end
t.StartTime = GetSecs;
FlushEvents('keyDown');
ListenChar;

GlobalTimer = 0; % this timer keeps track of all the timing in the experiment. TOTAL timing.
TimeUpdate = t.StartTime; % what time is it now?
%Present beginning of the begin fixation (with instructions)
Screen('FillOval', w, p.LeftRightColor(1,:), [CenterX-p.OuterFixationRadius-40 CenterY-p.OuterFixationRadius+25 CenterX+p.OuterFixationRadius-40 CenterY+p.OuterFixationRadius+25])
Screen('FillOval', w, p.LeftRightColor(2,:), [CenterX-p.OuterFixationRadius+20 CenterY-p.OuterFixationRadius+25 CenterX+p.OuterFixationRadius+20 CenterY+p.OuterFixationRadius+25])
Screen('DrawText',w, 'left', CenterX+p.OuterFixationRadius-38, CenterY-p.OuterFixationRadius*1.5+27, 20);
Screen('DrawText',w, 'right', CenterX+p.OuterFixationRadius+25, CenterY-p.OuterFixationRadius*1.5+27, 20);

Screen('Flip', w);
%TIMING!:
GlobalTimer = GlobalTimer + t.InstructionTime;
TimePassed = 0; % flush the time the previous event took
while (TimePassed<t.InstructionTime) % for as long as the begin fixation is on the screen...
    [resp, timestamp] = checkForResp([p.start,p.space],p.escape);
    if resp==-1; escaperesponse(OriginalCLUT); end 
    TimePassed = GetSecs-TimeUpdate; % and determine exactly how much time has passed since the start of the expt.
end
TimeUpdate = TimeUpdate + t.InstructionTime;
%Present the rest of the begin fixation
Screen('FillOval', w, p.FixColor, [CenterX-p.OuterFixationRadius CenterY-p.OuterFixationRadius CenterX+p.OuterFixationRadius CenterY+p.OuterFixationRadius])
Screen('Flip', w);
%TIMING!:
GlobalTimer = GlobalTimer + (t.BeginFixation-t.InstructionTime);
TimePassed = 0; % flush the time the previous event took
while (TimePassed<(t.BeginFixation-t.InstructionTime)) % for as long as the begin fixation is on the screen...
    [resp, timestamp] = checkForResp([p.start,p.space],p.escape);
    if resp==-1; escaperesponse(OriginalCLUT); end 
    TimePassed = GetSecs-TimeUpdate; % and determine exactly how much time has passed since the start of the expt.
end
TimeUpdate = TimeUpdate + (t.BeginFixation-t.InstructionTime);
ListenChar(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% A TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:p.NumTrials
% for n = 1:4
    
    t.TrialStartTime(n) = GlobalTimer; % get the starttime of each single trial (relative to experiment start)
    TimeUpdate = t.StartTime + t.TrialStartTime(n);
    p.TrialNumGlobal = p.TrialNumGlobal+1;
  
    %% Cue - subject must press now
    Screen('FillOval', w, p.LeftRightColor(p.WhichFinger(n),:), [CenterX-p.OuterFixationRadius CenterY-p.OuterFixationRadius CenterX+p.OuterFixationRadius CenterY+p.OuterFixationRadius]);
    Screen('DrawingFinished', w);
    t.FlipTimes(n,1) = GetSecs;
    Screen('Flip', w);
    
    %TIMING!:
    GlobalTimer = GlobalTimer + t.CueTime;
    CueTimePassed = 0; % flush time passed.
    respYet = 0;
    while (CueTimePassed<t.CueTime) % as long as the stimulus is on the screen...
        CueTimePassed = (GetSecs-TimeUpdate); % and determine exactly how much time has passed since the start of the expt.
        [resp, timestamp] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
        
        %scanner buttons are: b y g r (form left-to-right)
        if ~respYet && ~isempty(find(resp==p.keys,1))
            % response! 
            respYet=1;
            data.RespTimeRaw(n) = timestamp;
            data.RespTimeFromOnset(n) = timestamp - t.FlipTimes(n,1);
            data.Response(n) = find(resp==p.keys);    
            data.RespOnTime(n) = 1;
            if (find(resp==p.keys)==1 && p.WhichFinger(n)==1) || (find(resp==p.keys)==4 && p.WhichFinger(n)==2)
                data.Correct(n)=1;
            else
                data.WrongButton(n)=1;
            end              
        end
    end
    TimeUpdate = TimeUpdate + t.CueTime; % update Matlab on what time it is.
       
    %% ITI
    Screen('FillOval', w, p.FixColor, [CenterX-p.OuterFixationRadius CenterY-p.OuterFixationRadius CenterX+p.OuterFixationRadius CenterY+p.OuterFixationRadius]);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    %TIMING!:
    GlobalTimer = GlobalTimer + t.ITIs(n);
    itiTimePassed = 0; %Flush time passed.
    while (itiTimePassed<t.ITIs(n)) %As long as the stimulus is on the screen...
        itiTimePassed = (GetSecs-TimeUpdate); %And determine exactly how much time has passed since the start of the expt.
%         [keyIsDown, secs, keyCode] = KbCheck(-1);
        [resp, timestamp] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
        % checking to see if they respond during the ITI...this will count
        % as a "late" response for the previous trial
        if ~isempty(find(resp==p.keys,1))
            % response! 
            data.RespLate(n)=1;                      
        end
    end
    TimeUpdate = TimeUpdate + t.ITIs(n); %Update Matlab on what time it is.
  
end % end trial loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End fixation:
Screen('FillOval', w, p.FixColor, [CenterX-p.OuterFixationRadius CenterY-p.OuterFixationRadius CenterX+p.OuterFixationRadius CenterY+p.OuterFixationRadius]);
Screen('Flip', w);
GlobalTimer = GlobalTimer + t.EndFixation;
closingtime = 0; resp = 0;
while closingtime < t.EndFixation
    closingtime = GetSecs-TimeUpdate;
end
ListenChar(1); % unsuppressed keyboard mode
FlushEvents('keyDown'); % first discard all characters from the Event Manager queue
t.EndTime = GetSecs; % get endtime of the experiment in seconds
t.TotalExpTime = (t.EndTime-t.StartTime); % gets the duration of the run in secs.
t.GlobalTimer = GlobalTimer; % gets the duration of the run in secs.


%----------------------------------------------------------------------
%WINDOW CLEANUP--------------------------------------------------------
%----------------------------------------------------------------------
%This closes all visible and invisible screens and puts the mouse cursor
%back on the screen
Screen('CloseAll');
if exist('OriginalCLUT','var')
    if exist('ScreenNr','var')
        Screen('LoadCLUT', ScreenNr, OriginalCLUT);
    else
        Screen('LoadCLUT', 0, OriginalCLUT);
    end
end
clear screen
ShowCursor;



%----------------------------------------------------------------------
%PERFORMANCE CHECK-----------------------------------------------------
%----------------------------------------------------------------------

data.accuracy = mean(data.Correct);
nWrongTime = sum(data.RespLate);
nMissed = sum(~data.RespOnTime);
nWrongButton = sum(data.WrongButton);

fprintf('Overall accuracy: %.2f percent\n',data.accuracy*100);
fprintf('Number of trials where button pressed during ITI: %d\n',nWrongTime)
fprintf('Number of trials where failed to press a button: %d\n',nMissed)
fprintf('Number of trials where pressed wrong button: %d\n',nWrongButton)

%----------------------------------------------------------------------
%SAVE OUT THE DATA-----------------------------------------------------
%----------------------------------------------------------------------
cd(datadir);
%First I make a list of variables to save:
TheData(runnumber).t = t;
TheData(runnumber).p = p;
TheData(runnumber).data  = data;
eval(['save(''OriSpin_S', p.SubNum, '_DigitLocalizer.mat'', ''TheData'', ''-v7.3'')'])
cd(expdir)




