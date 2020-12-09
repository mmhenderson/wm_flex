% Get the timing of all events in all runs
% saves out a file that is used by MakeSampleFile

% MMH 9/5/18

clear
close all

% set inputs
FS_sbj = 'CP';
subnum = 7;
substr = sprintf('S%02d',subnum);

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
experiment_path = mypath(1:filesepinds(end-nDirsUp+1));

out_path = fullfile(experiment_path,'Samples');
beh_path = fullfile(experiment_path,'DataBehavior');

%% Timing information

TRdur = .8;
TRditched = 16; % 16 TRs were removed
not_recorded = TRdur * TRditched;


%% MAIN TASK %%%
currentRun = 0; nTRs_main = 583 - TRditched; 
% these are all going to be arrays [nTRs*nRunsTotal x 1], initialize them
% now.
RunLabels=[]; EventLabels = []; TrialLabels = []; RTLabels = [];
TargPos = []; BoundPos = []; CondLabels = []; CorrectResp = [];
RespActual = []; RandBoundPos = [];

trialnumglobal = 0;

nSess = 2;
nRuns = 10;

cols_odd = {[14.5922, 79.0745, 150.6510],[245.0980, 36.8980, 17.6039]};
cols_even = {[245.0980, 36.8980, 17.6039],[14.5922, 79.0745, 150.6510]};

bound_cols_odd_sess = {[153.5, 153.5, 153.5], [102.5, 102.5, 102.5]};
bound_cols_even_sess = {[102.5, 102.5, 102.5], [153.5, 153.5, 153.5]};

% get data file from each session (they are not concatenated here)
main_files = dir(fullfile(beh_path,substr,'Session*','*MainTask.mat'));
main_files = main_files(~contains({main_files.name},'TRAINING'));
assert(length(main_files)==nSess, 'should have exactly nSess sessions of main task data')

for sess = 1:nSess

%     load([beh_path 'S' char(subnum) filesep 'Session' num2str(sess) filesep main_files(sess).name])
    load(fullfile(main_files(sess).folder, main_files(sess).name));
  
    
    if sess==1
        nTrialsEach = TheData(1).p.NumTrials;
    end

    assert(length(TheData)==nRuns, 'should have %d runs per sess', nRuns)
   
    for rr = 1:nRuns

        currentRun = currentRun + 1;

        assert(contains(TheData(rr).p.Name, FS_sbj))
        % check some basic presentation parameters - colors of the
        % conditions and response sides. make sure all counterbalancing
        % was done correctly according to plan. 
        if mod(subnum,2)
            assert(all(TheData(rr).p.ConditionCueColors{1}==cols_odd{1}) &&...
                all(TheData(rr).p.ConditionCueColors{2}==cols_odd{2}));
        else
            assert(all(TheData(rr).p.ConditionCueColors{1}==cols_even{1}) &&...
                all(TheData(rr).p.ConditionCueColors{2}==cols_even{2}));
        end
        
        if mod(sess,2)
            assert(all(TheData(rr).p.BoundSideColors{1}==bound_cols_odd_sess{1}) &&...
                all(TheData(rr).p.BoundSideColors{2}==bound_cols_odd_sess{2}));
        else
            assert(all(TheData(rr).p.BoundSideColors{1}==bound_cols_even_sess{1}) &&...
                all(TheData(rr).p.BoundSideColors{2}==bound_cols_even_sess{2}));
        end
        
        % Find event times - there are four types of "events" here,
        % the target onset/offset and disk onsets/offsets. Make a couple of
        % lists nEvents long - describe the time (s) of the event, the type of
        % event, and the trial it corresponds to. Initialize these with
        % zeros.
        Event_onset = 0;  Event_type = 0; 
        Event_trial_global = 0; Event_trial_local = 0; % local is trial number in this run, global is over all sess.
   
        t = TheData(rr).t;
        p = TheData(rr).p;
        
        % list the dur of all events that happen in the trial...
        event_dur_list = [0, t.PreTargetCueTime, t.TargetTime, t.PreCueDelayLength, ...
            t.CueLength, t.BoundPreviewLength,t.PostCueDelayLength, t.BoundFinalLength];
        cum_event_dur_list = cumsum(event_dur_list);
        
        for n = 1:TheData(rr).p.NumTrials

            trialnumglobal = trialnumglobal+1;

            % list times for all events on this trial (seconds relative to
            % start of the nifti file)
            event_times_this_trial = cum_event_dur_list + t.TrialStartTime(n) - not_recorded;
            Event_onset = [Event_onset, event_times_this_trial];
            % pre-targ cue, targ, delay1, cue, bound-preview, delay2,
            % bound-actual, start ITI
            Event_type = [Event_type, [0.2, 1, 0, 0.3, 2, 0, 3, 0]];

            % mark the global and local trial number of these events
            Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
            Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];


        end
        
        % Find onset times for all TRs
        TR_onset = 0:.8:nTRs_main*0.8;
        TR_onset = TR_onset(1:end-1);

        % now convert sec to TRs
        % making a list of [nTRs x 1] for what happens on each TR
        triallabslocal_byTR = zeros(length(TR_onset),1);
        eventlabs_byTR = zeros(length(TR_onset),1);
        triallabsglobal_byTR = zeros(length(TR_onset),1);

        for i = 1:length(TR_onset)

            % on this TR - what type of of event is happening, and what trial is it a part of?
            % looking for the event label that was active during the first
            % half of this TR. If an event occurs in the second half of the
            % TR, then it gets marked in the next TR.
            middleofTR = TR_onset(i)+(TRdur/2);
            myind = find(middleofTR>Event_onset, 1, 'last');
            
            eventlabs_byTR(i) = Event_type(myind);
            triallabslocal_byTR(i) = Event_trial_local(myind);
            triallabsglobal_byTR(i) = Event_trial_global(myind);

        end

        % make sure that there is a 1 marking the start of every target
        % grating. This is necessary
        % because sometimes the events might not land in the right place in a
        % TR, because they are <0.8 s duration.
        for n = 1:nTrialsEach           
            eventlabs_byTR(find(triallabslocal_byTR==n & eventlabs_byTR~=0.2,1))= 1;
        end

        % now i'm appending to my long list of labels which will be saved. 
        numZeros = sum(triallabslocal_byTR==0); % any TRs that fell before first trial onset

        RunLabels = [RunLabels;repmat(currentRun,length(TR_onset),1)];
        TrialLabels = [TrialLabels; triallabsglobal_byTR];
        EventLabels = [EventLabels; eventlabs_byTR];

        % these things changed from trial to trial - indexing into these
        % arrays from my saved data to assign a value to every TR.
        TargPos = [TargPos; nan(numZeros,1); p.TargPosition(triallabslocal_byTR(numZeros+1:end))];
        BoundPos = [BoundPos;  nan(numZeros,1); p.BoundPosition(triallabslocal_byTR(numZeros+1:end))];
        CorrectResp = [CorrectResp;  nan(numZeros,1); p.CorrectResponse(triallabslocal_byTR(numZeros+1:end))];
        RespActual = [RespActual;  nan(numZeros,1); TheData(rr).data.Response(triallabslocal_byTR(numZeros+1:end))];
        RTLabels = [RTLabels; nan(numZeros,1); t.respTimeFromStart(triallabslocal_byTR(numZeros+1:end))];
        CondLabels = [CondLabels; nan(numZeros,1); p.Task(triallabslocal_byTR(numZeros+1:end))];
        if subnum~=1
            RandBoundPos = [RandBoundPos;  nan(numZeros,1); p.RandBoundPosition(triallabslocal_byTR(numZeros+1:end))];
        else
            RandBoundPos = [RandBoundPos; nan(length(triallabslocal_byTR),1)];
        end
    end %run loop

end

main.RunLabels = RunLabels;
main.TrialLabels = TrialLabels;
main.TargPos = TargPos;
main.BoundPos = BoundPos;
main.RespActual = RespActual;
main.CorrectResp = CorrectResp;
main.EventLabels = EventLabels;
main.CondLabels = CondLabels;
main.RTLabels = RTLabels;
main.RandBoundPos = RandBoundPos;

assert(numel(main.EventLabels)==nTRs_main*numel(unique(main.RunLabels)), 'wrong number of total TRs!')
assert(sum(main.EventLabels==1)==nRuns*nSess*nTrialsEach, 'make sure each trial is marked once')

%% SWM LOC TASK
if subnum~=1
    currentRun = 0; nTRs_swm = 508 - TRditched; 
    % these are all going to be arrays [nTRs*nRunsTotal x 1], initialize them
    % now.
    RunLabels=[]; EventLabels = []; TrialLabels = []; 
    TargPos = [];  RespActual = [];

    trialnumglobal = 0;

    nSess = 1;
    nRuns = 10;

    % get data file from each session (they are not concatenated here)
    swm_files = dir(fullfile(beh_path,substr,'Session*','*SWMLoc.mat'));
    swm_files = swm_files(~contains({swm_files.name},'TRAINING'));
    assert(length(swm_files)==nSess, 'should have exactly nSess sessions of swm task data')

    for sess = 1:nSess

        load(fullfile(swm_files(sess).folder, swm_files(sess).name));

        if sess==1
            nTrialsEach = TheData(1).p.NumTrials;
        end

        assert(length(TheData)==nRuns, 'should have %d runs per sess', nRuns)

        for rr = 1:nRuns
            
            currentRun = currentRun + 1;
            assert(contains(TheData(rr).p.Name, FS_sbj))
            % Find event times - there are four types of "events" here,
            % the target onset/offset and dial onset/offset. Make a couple of
            % lists nEvents long - describe the time (s) of the event, the type of
            % event, and the trial it corresponds to. Initialize these with
            % zeros.
            Event_onset = 0;  Event_type = 0; 
            Event_trial_global = 0; Event_trial_local = 0; % local is trial number in this run, global is over all sess.

            t = TheData(rr).t;
            p = TheData(rr).p;

            % list the dur of all events that happen in the trial...
            event_dur_list = [0, t.PreTargetCueTime, t.TargetTime, t.DelayLength, t.MaxRespTime];
            cum_event_dur_list = cumsum(event_dur_list);

            for n = 1:TheData(rr).p.NumTrials

                trialnumglobal = trialnumglobal+1;

                % list times for all events on this trial (seconds relative to
                % start of the nifti file)
                event_times_this_trial = cum_event_dur_list + t.TrialStartTime(n) - not_recorded;
                Event_onset = [Event_onset, event_times_this_trial];
                % pre-targ cue, targ, delay, resp period, start ITI
                Event_type = [Event_type, [0.2, 1, 0, 3, 0]];

                % mark the global and local trial number of these events
                Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
                Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];


            end

            % Find onset times for all TRs
            TR_onset = 0:.8:nTRs_swm*0.8;
            TR_onset = TR_onset(1:end-1);

            % now convert sec to TRs
            % making a list of [nTRs x 1] for what happens on each TR
            triallabslocal_byTR = zeros(length(TR_onset),1);
            eventlabs_byTR = zeros(length(TR_onset),1);
            triallabsglobal_byTR = zeros(length(TR_onset),1);

            for i = 1:length(TR_onset)

                % on this TR - what type of of event is happening, and what trial is it a part of?
                % looking for the event label that was active during the first
                % half of this TR. If an event occurs in the second half of the
                % TR, then it gets marked in the next TR.
                middleofTR = TR_onset(i)+(TRdur/2);
                myind = find(middleofTR>Event_onset, 1, 'last');

                eventlabs_byTR(i) = Event_type(myind);
                triallabslocal_byTR(i) = Event_trial_local(myind);
                triallabsglobal_byTR(i) = Event_trial_global(myind);

            end

            % make sure that there is a 1 marking the start of every target
            % grating. This is necessary
            % because sometimes the events might not land in the right place in a
            % TR.
            for n = 1:nTrialsEach           
                eventlabs_byTR(find(triallabslocal_byTR==n & eventlabs_byTR~=0.2,1))= 1;
            end

            % now i'm appending to my long list of labels which will be saved. 
            numZeros = sum(triallabslocal_byTR==0); % any TRs that fell before first trial onset

            RunLabels = [RunLabels;repmat(currentRun,length(TR_onset),1)];
            TrialLabels = [TrialLabels; triallabsglobal_byTR];
            EventLabels = [EventLabels; eventlabs_byTR];

            % these things changed from trial to trial - indexing into these
            % arrays from my saved data to assign a value to every TR.
            TargPos = [TargPos; nan(numZeros,1); p.TargPosition(triallabslocal_byTR(numZeros+1:end))];       
            RespActual = [RespActual;  nan(numZeros,1); TheData(rr).data.Response(triallabslocal_byTR(numZeros+1:end))];

        end %run loop

    end

    swm.RunLabels = RunLabels;
    swm.TrialLabels = TrialLabels;
    swm.EventLabels = EventLabels;

    swm.TargPos = TargPos;
    swm.RespActual = RespActual;

    assert(numel(swm.EventLabels)==nTRs_swm*numel(unique(swm.RunLabels)), 'wrong number of total TRs!')
    assert(sum(swm.EventLabels==1)==nRuns*nSess*nTrialsEach, 'make sure each trial is marked once')
else
    swm = [];
end

%% DWM LOC TASK
if subnum~=1
    currentRun = 0; nTRs_dwm = 452 - TRditched; 
    % these are all going to be arrays [nTRs*nRunsTotal x 1], initialize them
    % now.
    RunLabels=[]; EventLabels = []; TrialLabels = []; 
    ExpDigit = [];  ActDigit = []; RT = [];
    RespEarly = [];

    trialnumglobal = 0;

    finger_cols_1 = [1.5922 79.0745 155; 200 0 226.6; 80 80 80];
    finger_cols_0 = [200 0 226.6; 1.5922 79.0745 155; 80 80 80];

    
    % get data file from each session (they are not concatenated here)
    dwm_files = dir(fullfile(beh_path,substr,'Session*','*DWMLoc.mat'));
    dwm_files = dwm_files(~contains({dwm_files.name},'TRAINING'));
    fprintf('found %d sessions of DWM loc data\n',length(dwm_files));

    for sess = 1:length(dwm_files)

        load(fullfile(dwm_files(sess).folder, dwm_files(sess).name));

        if sess==1
            nTrialsEach = TheData(1).p.NumTrials;
        end

        for rr = 1:length(TheData)
            
            currentRun = currentRun + 1;
            
            assert(contains(TheData(rr).p.Name, FS_sbj))
            % check the response mapping that was used
            if (mod(subnum,2) && mod(rr,2)) || (~mod(subnum,2) && ~mod(rr,2))
                % odd sub, odd run, or even sub, even run
                assert(all(TheData(rr).p.LeftRightColor(:)==finger_cols_1(:)))
            else
                % even sub, odd run, or odd sub, even run
                assert(all(TheData(rr).p.LeftRightColor(:)==finger_cols_0(:)))
            end
            
            % Find event times - there are four types of "events" here,
            % the target onset/offset and dial onset/offset. Make a couple of
            % lists nEvents long - describe the time (s) of the event, the type of
            % event, and the trial it corresponds to. Initialize these with
            % zeros.
            Event_onset = 0;  Event_type = 0; 
            Event_trial_global = 0; Event_trial_local = 0; % local is trial number in this run, global is over all sess.

            t = TheData(rr).t;
            p = TheData(rr).p;

            % list the dur of all events that happen in the trial...
            event_dur_list = [0, t.TargetTime, t.DelayLength, t.MaxRespTime];
            cum_event_dur_list = cumsum(event_dur_list);

            for n = 1:TheData(rr).p.NumTrials

                trialnumglobal = trialnumglobal+1;

                % list times for all events on this trial (seconds relative to
                % start of the nifti file)
                event_times_this_trial = cum_event_dur_list + t.TrialStartTime(n) - not_recorded;
                Event_onset = [Event_onset, event_times_this_trial];
                % targ, delay, resp period, start ITI
                Event_type = [Event_type, [1, 0, 3, 0]];

                % mark the global and local trial number of these events
                Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
                Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];


            end

            % Find onset times for all TRs
            TR_onset = 0:.8:nTRs_dwm*0.8;
            TR_onset = TR_onset(1:end-1);

            % now convert sec to TRs
            % making a list of [nTRs x 1] for what happens on each TR
            triallabslocal_byTR = zeros(length(TR_onset),1);
            eventlabs_byTR = zeros(length(TR_onset),1);
            triallabsglobal_byTR = zeros(length(TR_onset),1);

            for i = 1:length(TR_onset)

                % on this TR - what type of of event is happening, and what trial is it a part of?
                % looking for the event label that was active during the first
                % half of this TR. If an event occurs in the second half of the
                % TR, then it gets marked in the next TR.
                middleofTR = TR_onset(i)+(TRdur/2);
                myind = find(middleofTR>Event_onset, 1, 'last');

                eventlabs_byTR(i) = Event_type(myind);
                triallabslocal_byTR(i) = Event_trial_local(myind);
                triallabsglobal_byTR(i) = Event_trial_global(myind);

            end

            % make sure that there is a 1 marking the start of every trial.
            for n = 1:nTrialsEach           
                assert(eventlabs_byTR(find(triallabslocal_byTR==n,1))==1);
    %             eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
            end

            % now i'm appending to my long list of labels which will be saved. 
            numZeros = sum(triallabslocal_byTR==0); % any TRs that fell before first trial onset

            RunLabels = [RunLabels;repmat(currentRun,length(TR_onset),1)];
            TrialLabels = [TrialLabels; triallabsglobal_byTR];
            EventLabels = [EventLabels; eventlabs_byTR];

            % these things changed from trial to trial - indexing into these
            % arrays from my saved data to assign a value to every TR.
            ExpDigit = [ExpDigit; nan(numZeros,1); p.WhichFinger(triallabslocal_byTR(numZeros+1:end))];       
            ActDigit = [ActDigit;  nan(numZeros,1); TheData(rr).data.Response(triallabslocal_byTR(numZeros+1:end))];
            RespEarly = [RespEarly; nan(numZeros,1); TheData(rr).data.RespTooEarly(triallabslocal_byTR(numZeros+1:end))];
            RT = [RT; nan(numZeros,1); TheData(rr).t.respTimeFromStart(triallabslocal_byTR(numZeros+1:end))];

        end %run loop

    end

    dwm.RunLabels = RunLabels;
    dwm.TrialLabels = TrialLabels;
    dwm.EventLabels = EventLabels;

    dwm.ExpDigit = ExpDigit;
    dwm.ActDigit = ActDigit;
    dwm.RespEarly = RespEarly;
    dwm.RT = RT;

    fprintf('found %d total runs of DWM task\n',numel(unique(dwm.RunLabels)));

    assert(numel(dwm.EventLabels)==nTRs_dwm*numel(unique(dwm.RunLabels)), 'wrong number of total TRs!')
else
    dwm = [];
end
%% SPATIAL POSITION LOCALIZER 
nTRs_spatloc = 391 - TRditched;

currentRun = 0; 
trialnumglobal = 0;

RunLabels = []; 
TrialLabels = [];
EventLabels = []; 
PosLabels = [];

% get data file from each session (they are not concatenated here)
loc_files = dir(fullfile(beh_path, substr,'Session*', '*sIEM_1D*'));

% looping over all runs across all sessions for this subject.
for rr = 1:length(loc_files)
    currentRun = currentRun + 1;
    
    load(fullfile(loc_files(rr).folder, loc_files(rr).name));

    if rr==1
        nTrialsEach = p.nTrials;
    end
    assert(contains(p.subName, FS_sbj))
    
    % where are the centers of wedges 1-24? 
    wedge_centers = 7.5:15:360; % these are in deg, clockwise from vertical.
    % where are the centers of wedges on each randomized trial?
    center_pos_each_trial = wedge_centers(p.trial_wedge)';
    
    % note that the 12.8 seconds have already been trimmed off here, this
    % is from the time right after the 12.8 seconds have been completed.
    Trial_onset = p.frameTime(p.trial_frame);
    % these are super exact times - if you rounded them they'd be about 3
    % sec apart. 
    
    % loop over trials and save a list of information about the events on
    % each trial.
    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;
    for n = 1:nTrialsEach

        trialnumglobal = trialnumglobal+1;

        % list times for all events per trial
        % here there is only 1 event per trial, the stim going on.
        event_times_this_trial = [Trial_onset(n)];
        Event_onset = [Event_onset, event_times_this_trial];
        Event_type = [Event_type,  1]; % on
        Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
        Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];

    end
    % also concatenating on an offset, because there were no gaps between
    % stims but there's a blank period at the end. 
    event_times_this_trial = Trial_onset(n)+p.stim_dur;
    Event_onset = [Event_onset, event_times_this_trial];
    Event_type = [Event_type,  0]; % off
    Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
    Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];

    
    % Find onset times for all my TRs
    TR_onset = 0:TRdur:nTRs_spatloc*TRdur;
    TR_onset = TR_onset(1:end-1);

    % Now we want to go through and map every event to the TR when it
    % happened. We'll save out a list nTRs long, describing which event and
    % trial number is occuring at that TR.
    eventlabs_byTR = nan(1,length(TR_onset));
    triallabslocal_byTR = nan(1,length(TR_onset));
    triallabsglobal_byTR = nan(1,length(TR_onset));

    for i = 1:length(TR_onset)

        middleofTR = TR_onset(i)+(TRdur/2);
        myind = find(middleofTR>Event_onset, 1, 'last');

        % on this TR - what type of of event is happening, and what global
        % trial number and run number is it a part of?
        eventlabs_byTR(i) = Event_type(myind);
        triallabslocal_byTR(i) = Event_trial_local(myind);
        triallabsglobal_byTR(i) = Event_trial_global(myind);
    end

     % make sure that the first TR of every trial is a 1
    for n = 1:nTrialsEach           
        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
    end

    numZeros = sum(triallabslocal_byTR==0);

    RunLabels = [RunLabels, repmat(currentRun,1,length(TR_onset))];
    TrialLabels = [TrialLabels, triallabsglobal_byTR];
    EventLabels = [EventLabels, eventlabs_byTR];   
    PosLabels = [PosLabels;nan(numZeros,1);center_pos_each_trial(triallabslocal_byTR(numZeros+1:end))];
    
end

spatLoc.RunLabels = RunLabels';
spatLoc.TrialLabels = TrialLabels';
spatLoc.PosLabels = PosLabels;
spatLoc.EventLabels = EventLabels';

assert(numel(spatLoc.EventLabels)==nTRs_spatloc*numel(unique(spatLoc.RunLabels)), 'wrong number of total TRs!')

%% DIGIT LOCALIZER %%%
nRuns = 0; nTRs_digloc = 399 - 16;
RunLabels = []; ExpDigit = []; ActDigit = [];
EventLabels = []; TrialLabels = []; RT = [];

% get data file from all session (they are concatenated)
loc_file = dir(fullfile(beh_path,substr,'Session*','*DigitLocalizer.mat'));
load(fullfile(loc_file.folder, loc_file.name))

finger_cols_1 = [1.5922 79.0745 155; 200 0 226.6; 80 80 80];
finger_cols_0 = [200 0 226.6; 1.5922 79.0745 155; 80 80 80];

if strcmp(substr,'S01')
    % for this subject, there was a recon problem with one of the digit
    % localizer runs, so we're not using that pair of runs. 
    runs2use = [1,2,3,4,7,8];  % in the behavior file, which runs should we USE?
    TheData = TheData(runs2use);
end
if strcmp(substr,'S06')
    % fix run 1
    % this is a run where subject reported pressing 2 instead of 4, but
    % was otherwise good at task
    TheData(1).data.Response(TheData(1).data.Response==2) = 4;
end
nTrialsEach = TheData(1).p.NumTrials;

trialnumglobal = 0;

% loop over all runs collected this session
for rr = 1:length(TheData)

    nRuns = nRuns + 1;

    % Find stimulus onset times
    Trial_onset = TheData(rr).t.TrialStartTime-not_recorded;
    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;
    
    assert(contains(TheData(rr).p.Name, FS_sbj))
    % check the response mapping that was used
    if rr==1
        last_cols = TheData(rr).p.LeftRightColor;
    else        
        assert(all(TheData(rr).p.LeftRightColor(1,:)==last_cols(2,:)) &&...
            all(TheData(rr).p.LeftRightColor(2,:)==last_cols(1,:)))
        last_cols = TheData(rr).p.LeftRightColor;
    end
   
    for n = 1:nTrialsEach
        
        trialnumglobal = trialnumglobal+1;
       
        % list times for 5 types of events per trial
        event_times_this_trial = [Trial_onset(n),Trial_onset(n)+TheData(rr).t.CueTime]; % cue onset
            
        Event_onset = [Event_onset, event_times_this_trial];
        Event_type = [Event_type,  1 0]; %orient multiplied by 1000 means it is a dial orientation. 
        Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,length(event_times_this_trial))];  
        Event_trial_local = [Event_trial_local, repmat(n,1,length(event_times_this_trial))];
        
    end
    
    % Find onset times for all my TRs
    TR_onset = 0:.8:nTRs_digloc*0.8;
    TR_onset = TR_onset(1:end-1);

    % Now we want to go through and map every event to the TR when it
    % happened. We'll save out a list nTRs long, describing which event and
    % trial number is occuring at that TR.
    eventlabs_byTR = nan(1,length(TR_onset));
    triallabslocal_byTR = nan(1,length(TR_onset));
    triallabsglobal_byTR = nan(1,length(TR_onset));
    
    for i = 1:length(TR_onset)
        
        middleofTR = TR_onset(i)+(TRdur/2);
        myind = find(middleofTR>Event_onset, 1, 'last');
        
        % on this TR - what type of of event is happening, and what global
        % trial number and run number is it a part of?
        eventlabs_byTR(i) = Event_type(myind);

        triallabslocal_byTR(i) = Event_trial_local(myind);
        triallabsglobal_byTR(i) = Event_trial_global(myind);
    end

    % make sure that the first TR of every trial is a 1
    for n = 1:nTrialsEach           
        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
    end
    
    numZeros = sum(triallabslocal_byTR==0);
            
    RunLabels = [RunLabels, repmat(nRuns,1,length(TR_onset))];
    TrialLabels = [TrialLabels, triallabsglobal_byTR];
    EventLabels = [EventLabels, eventlabs_byTR];
    
    
     
    digit_exp =  TheData(rr).p.WhichFinger;
    digit_act = TheData(rr).data.Response;
    digit_act(digit_act==4) = 2;

    rts = TheData(rr).data.RespTimeFromOnset;
        
    % switch the non-presses to label 3 (withhold)
    digit_act(isnan(digit_act)) = 3;
   
    RT = [RT,nan(numZeros,1),rts(triallabslocal_byTR(numZeros+1:end))'];
    ExpDigit = [ExpDigit,nan(numZeros,1),digit_exp(triallabslocal_byTR(numZeros+1:end))'];
    ActDigit = [ActDigit,nan(numZeros,1),digit_act(triallabslocal_byTR(numZeros+1:end))'];
    
end %run loop

% end
digLoc.RunLabels = RunLabels';
digLoc.TrialLabels = TrialLabels';
digLoc.RT = RT';
digLoc.ExpDigit = ExpDigit';
digLoc.ActDigit = ActDigit';
digLoc.EventLabels = EventLabels';

if numel(digLoc.EventLabels)~=nTRs_digloc*numel(unique(digLoc.RunLabels))
    error('wrong number of total TRs!')
end




%% Save timing file
if ~exist(out_path, 'dir'), mkdir(out_path); end
filename = fullfile(out_path,['TimingFile_', substr]);
fprintf('saving file to %s\n',filename);
save(filename, 'main', 'swm','dwm','spatLoc','digLoc','-v7.3');


