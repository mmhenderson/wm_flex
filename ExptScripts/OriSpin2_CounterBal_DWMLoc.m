function trialInfo = OriSpin2_CounterBal_DWMLoc(subNum)
% make a counterbalanced sequence of blocks over single session of 
% OriSpin2 (digit WM) task

% this function makes one block at a time, and gets called before each block.

% trialInfo is a structure holding all info for running the task
% (no randomization has to happen outside this function)

    %% insert subject number here before running script 
    % clear 
    if nargin<1
        subNum = '99';
    end
    
    mfilename('fullpath');
    rng('default')
    trialInfo.SubNum = subNum; 
    trialInfo.dateStr = datestr(now,'yymmdd');

    trialInfo.rndseed = round(sum(100*clock));
    rng(trialInfo.rndseed);

    %% define some parameters
    
    nBlocksToMake = 1;  % how many blocks to split trials into?
    nTrialsBlock = 20;
    
    %% calculate more numbers
   
    trialInfo.nTrialsBlock = nTrialsBlock;
    trialInfo.nBlocksToMake = nBlocksToMake;
   
    %% which finger to press on each trial
    digit_order_full = zeros(nBlocksToMake, nTrialsBlock);
   
    for bb=1:nBlocksToMake

        digit_order = repmat([1;2], nTrialsBlock/2,1);
        digit_order = digit_order(randperm(length(digit_order)));

        digit_order_full(bb,:) = digit_order;
    end
    
    trialInfo.WhichFinger = digit_order_full;
    
    %% other info about timing

    t.ITIrange = [1,5];
    
    iti_list = linspace(t.ITIrange(1),t.ITIrange(2), nTrialsBlock);
    trialInfo.ITI = zeros(nBlocksToMake,nTrialsBlock);
   
    for bb = 1:nBlocksToMake
        iti_list_shuff = iti_list(randperm(nTrialsBlock));
        trialInfo.ITI(bb,:) = iti_list_shuff;
    end
    
     %% calculate the total length of each run (different for each b/c delay periods are different)
    % all this timing info will go straight into the task script, just to
    % make sure we get the right number of TRs

    t.TargetTime = 1;   
    t.DelayLength = 12;    
    t.MaxRespTime = 1;

    t.BeginFixation = 13;
    t.EndFixation = 8;

    % total run length
    t.trialLength = t.TargetTime + t.DelayLength + t.MaxRespTime + mean(t.ITIrange);
    t.totalTime = t.BeginFixation + t.trialLength*nTrialsBlock + t.EndFixation;

    t.TR = 0.800;  
    t.totalTRsEachRun = ceil(t.totalTime/t.TR);

    %% save in structure to load
    
    trialInfo.t = t;

end

