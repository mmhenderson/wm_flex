function sIEM_1D_mapping_task()
% Retinotopy stimulus presentation script. Includes 3 different modes:
% meridian mapping (passive fixation), rotating wedge (attend for local
% contrast change), expanding ring (passive fixation).
% Serences Lab @ UCSD
% SYP 4/5/19 -  adapted to show wedges along a single ring.
%               fixation contrast change detection task.
% KA 4/8/19 - made checkResp function within this script; deleted
%             vestigial/commented out code; changed random seed to current
%             best practices (was outdated); added "p.debug mode" to not
%             wait full pre-time for scanner and not calc positions;
%             improved timing performance (tweaked shuffling, moved wedge
%             mask calculation outside of the main trial loop).
% TIMING FOR 24 1D WEDGES: 3s * 24 stim * 4 repetition 
% + 12.8s padding (16TR * 0.8s) + 12s end (15TR * 0.8s) = 312.8s = 391 TRs

clear all;
% dbstop if error
%%
%--------------------------------------------------------------------------
% Basic settings and dialog box
%--------------------------------------------------------------------------
p.screenWidth = 21;    % cm [scanner = 24, lab = 40]
p.viewDist = 47;       % cm [ scanner = 47, lab = 52];

% seed the generator (include this as a gui field in case want to repeat and exact sequence)
rng('default')
rng('shuffle')
p.rng_settings = rng; % Save the random seed settings!!
s = p.rng_settings.Seed;

% default answers for param prompt
defAns = {'CN','','0.65',num2str(s),'0.8','1'};

%get subject info
prompt = {'Subject Name', 'Run Number', 'Fixation Contrast (0-1)','Random Seed','TR (sec)', 'Fullscreen?'};

box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);
if length(box)==length(defAns)                  % all fields must be entered
    p.subName = char(box{1});
    p.scanNum = str2double(box{2});
    p.fixContrast = str2double(box{3});
    p.checkColors = 2; % 1 = red/green, 2 = black/white, 3 = blue/yellow, 4 = pink/blue, 5 = red/yellow
    p.rndSeed = str2double(box{4});
    p.TR = str2num(box{5});                    % TR in seconds
    p.fullScreen = str2num(box{6});             % fullscreen for scanner (1=fullscreen)
    p.train = 0;
    p.debug = 0; % if debugging, don't wait 16 TR's at the beginning, and print some stuff out
    p.no_quad_item_dups = 1; % 1 = no quadrant repeats, 2 = no item repeats
    p.gammacorrect = 0; 
else    %if cancel button or not enough input, then just bail
    error('Not the right number of inputs to GUI!')
    return
end

if ~p.train
    p.fMRI = 1;
else
    p.fMRI = 0;
end

if p.gammacorrect
    expdir = pwd;
    GeneralUseScripts = '/home/serencesadmin/Documents/MATLAB/Sunyoung/RRGeneralUseScripts';
    addpath(genpath(GeneralUseScripts));
end

% Set up the file name. 
c = clock;
p.month = [c(2)];
p.day = [c(3)]; 
fileName = sprintf('Data/%s_sIEM_1D_%d%d_r%d',p.subName,p.month,p.day,p.scanNum);
if ~exist('Data','dir')
    mkdir('Data');
end
if ~p.debug % DON'T OVERWRITE FILE UNLES WE ARE DEBUGGING!!
    if exist(fileName)
        Screen('CloseAll');
        msgbox('File already exists!', 'modal')
        return;
    end
end

%%
%--------------------------------------------------------------------------
% Stimulus Parameters (Size, timing, color)
%--------------------------------------------------------------------------
% pre-make some color combos:
if p.checkColors == 1       % red green
    p.color1 = [1 0 0]; p.color2 = [0 1 0];
elseif p.checkColors == 2   % white and black
    p.color1 = [0 0 0]; p.color2 = [1 1 1];
elseif p.checkColors == 3   % blue yellow
    p.color1 = [0 0 1]; p.color2 = [1 1 0];
elseif p.checkColors == 4   % pink blue
    p.color1 = [0.5 0 0.5]; p.color2 = [0 0 0.5];
elseif p.checkColors == 5   % red yellow
    p.color1 = [0.5 0 0]; p.color2 = [1 1 0];
end

% Ideally we should update this to be exactly divisible by the refresh
% rate!!!!
p.tf = 4;           % Hz for FULL cycle (checkerboard reversal) so 4 Hz here for full cycle = 8 Hz for flashing checkerboard..)
p.fixSize = 10;     %fixation point diameter in pixels

% MAKE SURE SUBJECTS DON'T MOVE DURING THIS TIME
p.preDur = 16*p.TR; % junk before stimulus in seconds (throw away these) multi-band data REQUIRES first 16 TRs to be thrown out in order to properly reconstruct the data!
p.endDur = 15*p.TR; % blank at the end for hemodynamic delay

% PSY -- discrete mapping stimuli
p.nWedge = 24;
p.nRing = 1;
p.nRep = 4; % 4 repetitions of all stim in 1 run
p.nStim = p.nWedge*p.nRing*p.nRep; % 24 different stimuli * 3 repetition for each

% checkerboards should be same across run types.
p.nAng = p.nWedge*2;	% KA make number of sub-wedgees a multiple of the total # (So they all look the same!)
p.nRad = 15;	% number of sub-rings in a circle

% need to specify order randomization, indexing
p.meriThick = 1/(p.nWedge);
p.fixedRingRad = 7; % in deg
p.fixedRingThick = 2.5*2; % in deg
p.fixedRingIn = p.fixedRingRad-p.fixedRingThick/2;
p.fixedRingOut = p.fixedRingRad+p.fixedRingThick/2;

p.stim_dur = 3; % stim duration in sec
p.period = p.stim_dur*p.nStim; % each stimulus * 4 s of display?
p.innerRing = 0.03;
p.nCycles = 1;

p.startPhase = 0;
% index of all flips across the entire experiment, resets every stim cycle
flipIndex = repmat(1:p.period*p.tf*2, 1, p.nCycles);
stimIndex = ceil(flipIndex/(p.tf*2*p.stim_dur)); % shuffle later for each cycle

%%
%--------------------------------------------------------------------------
% RESPONSE PARAMS
%--------------------------------------------------------------------------
KbName('UnifyKeyNames');
p.escape = KbName('ESCAPE');
p.space = KbName('space');

if p.fMRI
    p.start = KbName('t');  % for scanner
    p.keys = [KbName('b'),KbName('1')]; % only response key, for detection of contrast change -- put in boht letters and numbers just in case
    p.wrong_keys = [KbName('y'),KbName('r'),KbName('g'),KbName('2'),KbName('3'),KbName('4')]; % check if they pressed the wrong key.... 
else
    p.start = p.space; % subject can press enter to start (make sure this is the same for windows)
    p.nextTrial = p.space;
    p.keys = KbName('b');  %'1!'); % only response key, for detection of contrast change
end
ListenChar(2); %make sure key presses aren't entered into this script while it's running.
%%
%--------------------------------------------------------------------------
% OPEN THE WINDOW, SCREEN PARAMS
%--------------------------------------------------------------------------
% figure out how many screens we have, and pick the last one in the list
s = max(Screen('Screens'));
Screen('Preference', 'SkipSyncTests', 1);
% compute middle gray
p.black = BlackIndex(s);
p.white = WhiteIndex(s);
p.gray=round((p.white+p.black)/2);
if round(p.gray)==p.white
    p.gray = p.black;
end

if p.fullScreen
    [wptr, ScreenSize] = Screen('OpenWindow', s, p.gray);
else
    % if we're dubugging open a 640x480 window that is a little bit down from the upper left
    % of the big screen
    %     [wptr, ScreenSize]=Screen('OpenWindow', s, p.gray, [50,50,1024+50,768+50]);
    %     [wptr, ScreenSize]=Screen('OpenWindow', s, p.gray, [50,50,50+50,50+50]);
    [wptr, ScreenSize]=Screen('OpenWindow', s, p.gray, [0 0 1024 768]);
end

%----------------------------------------------------------------------
%WINDOW SETUP & GAMMA CORRECTION---------------------------------------
%----------------------------------------------------------------------
if p.gammacorrect
    PsychJavaTrouble;

    OriginalCLUT = Screen('LoadClut', wptr);
    MyCLUT = zeros(256,3); MinLum = 0; MaxLum = 1;
    CalibrationFile = 'calib_07-Nov-2016.mat';
    [gamInverse,dacsize] = LoadCalibrationFileRR(CalibrationFile, expdir, GeneralUseScripts);
    LumSteps = linspace(MinLum, MaxLum, 256)';
    MyCLUT(:,:) = repmat(LumSteps, [1 3]);
    MyCLUT = round(map2map(MyCLUT, repmat(gamInverse(:,4),[1 3]))); %Now the screen output luminance per pixel is linear!
    Screen('LoadCLUT', wptr, MyCLUT);
    clear CalibrationFile gamInverse
end

% Screen timing information
p.hz = Screen('NominalFrameRate', wptr); % get refresh rate in Hz juuust in case!
p.ifi = Screen('GetFlipInterval', wptr); % get inter-flip interval (i.e. the refresh rate) -- Note, this is more precisely measured by PTB than "frame rate", so use this for timing sensitive things!
p.aspectRatio = ScreenSize(4)/ScreenSize(3);

% figure out pixel size and stuff
p.pxSize = p.screenWidth ./ ScreenSize(3);          % each pixel is this many cm wide
screenHeight_DVA = rad2deg(2*atan( (p.screenWidth*p.aspectRatio) ./ (2*p.viewDist)));
screenWidth_DVA = screenHeight_DVA / p.aspectRatio;

% Check GPU stuff
AssertOpenGL;
AssertGLSL;

% Enable alpha blending for typical drawing of masked textures:
Screen('BlendFunction', wptr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Create a special texture drawing shader for masked texture drawing:
glsl = MakeTextureDrawShader(wptr, 'SeparateAlphaChannel');

if p.fullScreen
    HideCursor;	% Hide the mouse cursor
    % set the priority up way high to discourage interruptions
    Priority(MaxPriority(wptr));
end

%%
%--------------------------------------------------------------------------
% MAKE STIMULI
%--------------------------------------------------------------------------
% fixation rectangle
fixRect = [(ScreenSize(3)-p.fixSize)/2, (ScreenSize(4)-p.fixSize)/2,...
    (ScreenSize(3)+p.fixSize)/2, (ScreenSize(4)+p.fixSize)/2 ];

% make cartesian matries
[x,y] = meshgrid( linspace(-1/p.aspectRatio,1/p.aspectRatio,ScreenSize(3)), ...
    linspace(-1,1,ScreenSize(4)));

% make polar angle matrices
r = sqrt(x.^2+y.^2);        % distance from center, in aspect ratio
a = atan2(y,x)/(2*pi)+.5;   % 0<=a<=1

% hm...let's compute DVA polar angle matrices
yd = linspace(-screenHeight_DVA/2,screenHeight_DVA/2,ScreenSize(4));
xd = linspace(-screenWidth_DVA/2,screenWidth_DVA/2,ScreenSize(3));
[xx,yy] = meshgrid(xd,yd);
rd = sqrt(xx.^2+yy.^2);
ad = atan2(yy,xx)/(2*pi)+.5;   % 0<=a<=1

% Duncan & Boynton 2003 cortical magnification factor
M_db = (0.065*rd + 0.054);          % this gives us cortical mm at this deg
norm_cmf = (-1*M_db) + max(abs(M_db(:))) + 1;

%make the checkerboard
radCheck = sign( sin( norm_cmf.*r.*p.nRad.*pi ) .* sin(p.nAng.*pi.*ad) ); % original

% Assign color to matrix for later
rgbcheck = zeros(ScreenSize(4)*ScreenSize(3), 3, 'uint8');
rgbcheck2 = zeros(ScreenSize(4)*ScreenSize(3), 3, 'uint8');
for col = 1:3
    indx = find(radCheck(:) == 1);	 % make these pixels one color
    rgbcheck(indx, col) = p.color1(col); % fill in the red/gren/blue plane, scale by desired red channel value
    % then make the other colored checkerboard (counterphase to the first)
    rgbcheck2(indx, col) = p.color2(col);
    % then do the same, but for the second color
    indx = find(radCheck(:) == -1);
    rgbcheck(indx, col) = p.color2(col);
    rgbcheck2(indx, col) = p.color1(col);
end

% reshape into a XxYx3
rgbcheck = reshape(rgbcheck, [ScreenSize(4), ScreenSize(3), 3]);
rgbcheck2 = reshape(rgbcheck2, [ScreenSize(4), ScreenSize(3), 3]);

% make and draw your two background checkerboards (these will alternate to
% make a flip...
textureIndex(1)=Screen('MakeTexture', wptr, (rgbcheck)*255);
textureIndex(2)=Screen('MakeTexture', wptr, (rgbcheck2)*255);

% initialize the masks - realmask will equal a grey screen + an alpha
% channel that will 'reveal' the checkerboard (hence the 2nd dim)
realmask = ones(ScreenSize(4), ScreenSize(3),2, 'uint8')*p.gray;	% middle gray
mask = zeros(ScreenSize(4), ScreenSize(3), 'uint8');                % fully transparent mask

% fixation brightness change
p.fix1 = 0;
p.fix2 = p.gray*p.fixContrast;
p.nFixChange = 20; % EVEN NUMBER!!! - # of changes for a run

% brightness change -- assign change points with at least # flips apart
cutflips = floor(length(flipIndex)/(p.nFixChange*2));
p.changepoint = randi(cutflips,[1 p.nFixChange]);
for i = 1:p.nFixChange
    p.changepoint(i) = p.changepoint(i)+(i-1)*cutflips*2+floor(cutflips/2);
end
p.min_change_time = min(diff(p.changepoint)); % minimum number of "frames" between each fiation point change

% Add the last timepoint to the "changepoint" vector so that we can score
% accuracy for the last change!!
p.check_changepoint_resp = [p.changepoint,numel(flipIndex)];
% Reshape this vector into 2 columns so that every other change is from light to dark, back to
% light again, etc.
p.changepoint2 = (reshape(p.changepoint,[2,p.nFixChange/2]))';
p.fixmat = ones(1,length(flipIndex));
for i = 1:size(p.changepoint2,1)
    p.fixmat((p.changepoint2(i,1)):(p.changepoint2(i,2))) = 2;
end

% % PSY: make semi-random presentation order
% find three segments and piece together
Screen('TextSize', wptr, 25);
DrawFormattedText(wptr,'Loading...This might take a minute.', 'center', ScreenSize(4)/2.5);
Screen('FillRect', wptr, [0 0 0], fixRect );
Screen('Flip', wptr);
%% DON'T ALLOW LOCATION REPEATS (QUADRANT OR EXACT)
justonerep = p.nWedge*p.nRing;

% Preallocate some matrices to make our shuffling a lot faster
saveOrder = NaN(justonerep,p.nRep); % made reshaping easier by putting repetitions in columns!
saveShuffled = NaN(justonerep,p.nRep);

% For 24 wedges, 1-6 = upper right, 7 - 12 = lower right, 13-18 = lower
% left, 19- 24 = upper left.
all_locs = 1:justonerep; all_locs = all_locs';
all_quads = [ones(1,justonerep/4),ones(1,justonerep/4)*2,ones(1,justonerep/4)*3,ones(1,justonerep/4)*4]; all_quads = all_quads';
shuffle_ind = 1:justonerep;
countO = 1;
while countO < p.nRep+1 % this code shuffles presentation of wedges so that we don't have any quadrant repeats (though do we not want any? this makes quadrants temporally anticorrelated?)
    perm_ind = Shuffle(shuffle_ind);
    shuffled_locs = all_locs(perm_ind);
    shuffled_quads = all_quads(perm_ind);
    if p.no_quad_item_dups == 1 % don't allow any back-to-back items in the same quadrant
        if sum(abs(diff(shuffled_quads))==0)==0 % ~length(find(diff(temp)==0)) % replaced this because it was slower (KA)
            if countO == 1
                saveOrder(:,countO) = shuffled_quads;
                saveShuffled(:,countO) = shuffled_locs;
                countO = countO+1;
            elseif sum(abs(saveOrder(:,countO-1)-shuffled_quads))~=0 && (saveOrder(end,countO-1) ~= shuffled_quads(1)) % made this a joint logical because matlab didn't like the .* with logicals
                saveOrder(:,countO) = shuffled_quads;
                saveShuffled(:,countO) = shuffled_locs;
                countO = countO+1;
            end
        end
    else % just don't allow exact position repeats
        if sum(abs(diff(shuffled_locs))==0)==0  % ~length(find(diff(temp)==0)) % replaced this because it was slower (KA)
            if countO == 1
                saveOrder(:,countO) = shuffled_quads;
                saveShuffled(:,countO) = shuffled_locs;
                countO = countO+1;
            elseif sum(abs(saveOrder(:,countO-1)-shuffled_locs))~=0 && (saveOrder(end,countO-1) ~= shuffled_locs(1)) % made this a joint logical because matlab didn't like the .* with logicals
                saveOrder(:,countO) = shuffled_quads;
                saveShuffled(:,countO) = shuffled_locs;
                countO = countO+1;
            end
        end
    end
end

p.quadOrder = saveOrder(:); % quadrant order randomized
p.randOrder = saveShuffled(:); % specific location order randomized
%% DRAW ALL THE APERTURE TEXTURES OUTSIDE THE MAIN LOOP! 1 texture per wedge stimulus.
for w = 1:p.nWedge
    this_stim = w;
    % mask stuff
    mask = mask.*0; % zero out the masks (faster than reallocating each time)
    realmask = realmask.*0;
    realmask = (realmask+1)*p.gray;  %turn it gray
    % basic circle with aperture
    mask(r<p.innerRing)=255;
    mask(r>max(y(:)))=255;
    
    % PSY
    p.angLo = p.meriThick*(this_stim-1)+0.25; % +0.25 to shift the wedges, align to vertical
    p.angHi = p.meriThick*this_stim+0.25;
    
    if p.angLo < 0
        mask(a > p.angHi & a<p.angLo+1) = 255;
    elseif p.angHi > 1 && p.angLo < 1
        mask(a>p.angHi-1 & a<p.angLo) = 255;
    elseif p.angHi > 1 && p.angLo >= 1
        mask(a<p.angLo-1) = 255;
        mask(a>p.angHi-1) = 255;
    else
        mask(a<p.angLo) = 255;
        mask(a>p.angHi) = 255;
    end
    
    % 1D mapping around single fixed ring -- no need for CMF!
    radLo = p.fixedRingIn;
    radHi = p.fixedRingOut;
    
    mask(rd < radLo )=255;
    mask(rd > radHi )=255;
    
    realmask(:,:,2) = mask;
    maskTexture{w} = Screen('MakeTexture', wptr, realmask, [], [], [], [], glsl); % make this have an alpha channel
    
end
%% SETUP SPACE FOR VARIABLES
% define some counters & such for trial loop (Do this before starting the
% scanner timing)
targCntr = 0;           % number of targets
nextFlip = 0;           % time that the next image should be drawn
p.frameTime = NaN(numel(flipIndex),1);
% p.phasePerFrame = NaN(numel(flipIndex),1);
p.keyResp = zeros(numel(flipIndex),1);  % did they press a key on this frame?
p.wrongKeyResp = zeros(numel(flipIndex),1);  % did they press the wrong key on this frame? Save it! We can always 

p.nTrials = p.nWedge*p.nRing*p.nRep;
p.trial = NaN(p.nTrials,1);
p.trial_frame = NaN(p.nTrials,1);
p.trial_wedge = NaN(p.nTrials,1);
p.trial_fix_color = NaN(p.nTrials,1);

%% BEGIN DISPLAY

Screen('TextSize', wptr, 25);
centerY = ScreenSize(4)/2; 
DrawFormattedText(wptr,'Press button when brightness of the fixation dot changes.', 'center', centerY - 100);
DrawFormattedText(wptr,'Please stay perfectly still, even when the screen is blank.', 'center', centerY - 50);
DrawFormattedText(wptr,'Scanner Ready! Press start key to continue.', 'center', centerY + 100);
if p.debug
    debugText = 'Warning! You are in debug mode!';
    DrawFormattedText(wptr,debugText, 'center', centerY - 100);
end
Screen('FillRect', wptr, [0 0 0], fixRect );
Screen('Flip', wptr);

% WAIT FOR SPECIFIC KEYPRESS TO START (e.g. scanner trigger)
% check for quit response
resp = 0;
while resp == 0
    [resp, ~] = checkForResp([p.start]);
end

FlushEvents;

% Draw a blank screen for the pre-task 16 TR's 
Screen('FillRect', wptr, [0 0 0], fixRect );
Screen('Flip', wptr);


if p.debug
    fprintf('trial\twedge\nflip\tfix\n');
end

% Wait a designated amount of time before starting!
if ~p.train
    % wait some amount of time before starting.
    if p.debug
        p.preDur = 1;
        WaitSecs(p.preDur);
    else
        WaitSecs(p.preDur);
    end
end

% Start timing RIGHT before we start the trial loop (don't move!)
tic
startTime = GetSecs;
%% TRIAL LOOP
for flipCount=1:numel(flipIndex)
    
    
    %% TARGET PRESENTATION
    this_stim = p.randOrder(stimIndex(flipCount)); % stimulus to present
    
    if p.fixmat(flipCount) == 1
        fixcol = p.fix1;
    elseif p.fixmat(flipCount) == 2
        fixcol = p.fix2;
    end
    % print out this trial if debugging and/or save some more variables!
    if flipCount == 1
        targCntr = targCntr+1;
        if p.debug % don't print stuff if we're not debugging
            fprintf('%d\t%d\t%d\t%d\n', targCntr, this_stim, flipCount, p.fixmat(flipCount));
        end
        % Save for easy analysis!
        p.trial(targCntr) = targCntr;
        p.trial_frame(targCntr) = flipCount;
        p.trial_wedge(targCntr) = this_stim;
        p.trial_fix_color(targCntr) = p.fixmat(flipCount);
    else
        if this_stim ~= p.randOrder(stimIndex(flipCount-1))
            targCntr = targCntr+1;
            if p.debug
                fprintf('%d\t%d\t%d\t%d\n', targCntr, this_stim, flipCount, p.fixmat(flipCount));
            end
            % Save for easy analysis!
            p.trial(targCntr) = targCntr;
            p.trial_frame(targCntr) = flipCount;
            p.trial_wedge(targCntr) = this_stim;
            p.trial_fix_color(targCntr) = p.fixmat(flipCount);
        end
    end
    
    
    %% STIMULUS MASK
    
    % now draw the background checkerboard, then the mask
    Screen('DrawTexture',wptr,textureIndex(mod(flipCount,2)+1));
    Screen('DrawTexture',wptr,maskTexture{this_stim});
    Screen('FillRect', wptr, fixcol, fixRect);
    Screen('DrawingFinished',wptr);
   
   
    
    %% RESPONSE FLIP THIS FRAME
    % wait to flip it
    resp = 0; 
    while toc < nextFlip %%% moved checkForResp in here to continuously check response while waiting for the next thing 
        
        % check response...
        [resp, ~] = checkForResp([p.keys, p.wrong_keys, p.escape]); % checks both buttons...
        if GetSecs==-1
            ListenChar(0);
            sca;
            save(fileName,'p');
            return
        end
        
        if resp && any(p.keys==resp) %if subject responds with a button press
            p.keyResp(flipCount) = 1;
        elseif resp && any(p.wrong_keys==resp)
            p.wrongKeyResp(flipCount) = 1;
        end
        
        if resp==p.escape
            ListenChar(0);
            Screen('CloseAll');
            save(fileName,'p');
            return
        end
        
        
    end
    
    % do the flip
    Screen('Flip',wptr);
    nextFlip = flipCount/(2*p.tf); %time at which the next flip should be performed
    
    p.frameTime(flipCount) = GetSecs - startTime;
    
end

%% END OF EXPERIMENT - TOTAL TIME
%wait for last frame to finish
while toc < nextFlip
end
endTime = GetSecs;

Screen('FillRect', wptr, [0 0 0], fixRect);
Screen('Flip', wptr);
WaitSecs(p.endDur); 

p.totalTaskTime = endTime - startTime; % just the task time.. .should be 72 stimuli x 4 sec each = 288
p.totalRunTime = p.totalTaskTime + p.preDur + p.endDur; % Time including the pre trial baseline!
p.totalTRs = p.totalRunTime / p.TR;

save(fileName,'p'); % save once in case we crash it doing something dumb while calculating acc

%% CALCULATE ACCURACY
% let's introduce a response deadline!

p.mean_frame_time = nanmean(diff(p.frameTime));
% minimum between stimuli with 100 changes = 1.48 seconds...
% Do 2 seconds or whichever was smaller.....
p.resp_deadline = min([ round2(2,p.mean_frame_time),  round2(min(p.changepoint)*p.mean_frame_time,p.mean_frame_time)]); % seconds to respond -- rounded to nearest "flikcer frame" where we collected response
p.resp_min = round2(.15,p.mean_frame_time); % False alarm type 1 (responded in <100 ms or so)

nFrames_max = p.resp_deadline / p.mean_frame_time;
nFrames_min = p.resp_min / p.mean_frame_time;

p.Hit = NaN(length(p.changepoint),1);
p.FA = NaN(length(p.changepoint),1);
p.MultipleResponses = NaN(length(p.changepoint),1);

for t = 1:length(p.changepoint)
    % this checks if they detected the target at ALL before the next
    % change.. but sometimes this is a fairly lonnnngg amount of time...
    % sometimes any seconds.
    resp_during_targ = sum(p.keyResp(p.check_changepoint_resp(t)+nFrames_min: p.check_changepoint_resp(t)+nFrames_max))~=0;
    resp_during_nontarg = sum(p.keyResp([p.check_changepoint_resp(t):p.check_changepoint_resp(t)+nFrames_min-1,p.check_changepoint_resp(t)+nFrames_max+1:p.check_changepoint_resp(t+1)]))~=0;
    
    p.Hit(t) = resp_during_targ;
    p.FA(t) = resp_during_nontarg;
    p.MultipleResponses(t) = resp_during_targ & resp_during_nontarg;
end
p.runHits = sum(p.Hit)/p.nFixChange;
p.runFA = sum(p.FA)/p.nFixChange;
p.runMultipleResponses = sum(p.MultipleResponses)/p.nFixChange;

% fprintf('Accuracy: %2.1f\n',p.runAcc);

%% END EXPERIMENT & CLEANUP
outputString = 'All done!';
accString = sprintf('Hit Rate: %2.1f\n',p.runHits);
faString = sprintf('False Alarm Rate: %2.1f\n',p.runFA);
multipleString = sprintf('Too many responses!: %2.1f times\n',p.runMultipleResponses);
wrongkeyString = sprintf('Wrong key pressed! %2.1f times\n',sum(p.wrongKeyResp));

Screen('TextSize', wptr, 20);
Screen('DrawText', wptr, outputString, [ScreenSize(3)/4], [ScreenSize(4)/2]);
Screen('DrawText', wptr, accString, [ScreenSize(3)/4], [ScreenSize(4)/2]+50);
Screen('DrawText', wptr, faString, [ScreenSize(3)/4], [ScreenSize(4)/2]+75);
Screen('DrawText', wptr, multipleString, [ScreenSize(3)/4], [ScreenSize(4)/2]+100);
Screen('DrawText', wptr, wrongkeyString, [ScreenSize(3)/4], [ScreenSize(4)/2]+125);

Screen('Flip', wptr);


save(fileName,'p');

KbWait; %wait for keypress
ListenChar(0);
ShowCursor; % show cursor again if we hid it (full screen mode)

% if gammacorrect
%     Screen('LoadCLUT', wptr, OriginalCLUT);
% end
Screen('CloseAll');

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Check for a response
%--------------------------------------------------------------------------
function [resp, ts] = checkForResp(possResp)
% queries the keyboard to see if a legit response was made
% returns the response and the timestamp
resp = 0;
ts = nan;
[keyIsDown, secs, keyCode] = KbCheck;

if sum(keyCode)>=1   % if at least one key was pressed
    keysPressed = find(keyCode);
    % in the case of multiple keypresses, just consider the first one
    if find(keysPressed(1)==possResp)
        resp = keysPressed(1);
        ts = secs;
    end
end
end
%--------------------------------------------------------------------------
% Round to nearest multiple of something
%--------------------------------------------------------------------------
function z = round2(x,y)
%% defensive programming
error(nargchk(2,2,nargin))
error(nargoutchk(0,1,nargout))
if numel(y)>1
    error('Y must be scalar')
end

%%
z = round(x/y)*y;
end

%--------------------------------------------------------------------------
% Loading gamma file
%--------------------------------------------------------------------------
function [gamInverse,dacsize] = LoadCalibrationFileRR(CalibrationFile, expdir, GeneralUseScripts)

%LoadCalibrationFile.m is a function that loads and returns gamInverse & dacsize 
%from the calibrationFile.
%This script is created by FT on 2000/05/10 and modified by JJ on 2007/06/28
%Then, RR made this version on 2010/10/10 so that it will work on a pc, 
%where paths sopmetimes contain spaces. 

err = 0;

eval(['cd ''' GeneralUseScripts '/Calibration''']);

eval(['load ' CalibrationFile ' gamInverse dacsize'],'err = 1'); %Load gamInverse & dacsize from CalibrationFile

if err
	disp(['Warning! Calibration File Not Found!']);
	query = input('Do you want to proceed using uncorrected gamma?  y or n(default) ??  ','s');
	if ~isempty(query) && query(1) == 'y'
		gamInverse = [0:255; 0:255; 0:255]'; %Use uncorrected gamma inverse function
		dacsize = 8; %Use dacsize 8
	else
		error('Exiting program...');
	end
end
eval(['cd ''' expdir '''']);
end
