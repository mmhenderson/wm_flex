function trialInfo = OriSpin2_CounterBal_SWMLoc(subNum,makePlots)
% make a counterbalanced sequence of blocks over single session of 
% OriSpin2 (spatial WM) task

% trialInfo is a structure holding all info for running the task
% (no randomization has to happen outside this function)

    %% insert subject number here before running script 
    % clear 
    if nargin<1
        subNum = '99';
        makePlots=0;
    end
    
    mfilename('fullpath');
    rng('default')
    trialInfo.SubNum = subNum; 
    trialInfo.dateStr = datestr(now,'yymmdd');

    trialInfo.rndseed = round(sum(100*clock));
    rng(trialInfo.rndseed);
    
    %% define some parameters
    nSess = 1;
    nBlocksTotal = 10;  % how many blocks to split trials into?
    nTrialsBlock = 20;
    nBinsTarg = nTrialsBlock; % How mamy bins do we divide the stim space into?
    nDegrees = 360;
    axis_spacing = nDegrees/(nTrialsBlock*nBlocksTotal);
    % define all the possible spatial positions that can be sampled
    orient_space_targ = axis_spacing:axis_spacing:nDegrees;
    % break this space into bins 
    binListTarg = reshape(orient_space_targ, ...
        length(orient_space_targ)/nBinsTarg, nBinsTarg);
    
    %% calculate more numbers
    nBlocksSess = nBlocksTotal/nSess;

    trialInfo.nSess = nSess;
    trialInfo.nTrialsBlock = nTrialsBlock;
    trialInfo.nBlocksSess = nBlocksSess;
   
    %% do all the counter-balancing here
    % initalize these things over all trials at once
    targ_position = zeros(nSess,nBlocksSess, nTrialsBlock);
    % decide where the continuous report dot will start on each trial (random)
    start_position = zeros(nSess, nBlocksSess, nTrialsBlock);
    
    for ss = 1:nSess
        % randomly assign one stim from each bin to each run.
        % shuffling within each bin (columns of this matrix).
        % afterward, each row of this matrix will define the stims to use on a
        % single run. 
        targ_pos_order_rand = binListTarg;
        for bin=1:nBinsTarg
            targ_pos_order_rand(:,bin) = targ_pos_order_rand(randperm(nBlocksTotal),bin);
        end
        
        % same thing for the response probe starting positions
        % sample once from each bin each run. 
        start_pos_order_rand = binListTarg;
        for bin=1:nBinsTarg
            start_pos_order_rand(:,bin) = start_pos_order_rand(randperm(nBlocksTotal),bin);
        end
        
        for bb = 1:nBlocksSess       
            % already know which stims to use on this block, but not yet
            % their order. Define a random order for the trials here
            targ_position(ss,bb,:) = targ_pos_order_rand(bb,randperm(nBinsTarg));
            % start positions and target positions are shuffled
            % independently, so the correspondence between them is random.
            start_position(ss,bb,:) = start_pos_order_rand(bb,randperm(nBinsTarg));           
        end
    end
    
    % these are the final variables we will reference during the experiment
    trialInfo.TargPosition = targ_position;
    trialInfo.ReportStartPosition = start_position;
    %% double check my counter-balancing
    % making sure that each run includes one target and one response
    % starting position from each bin.
    count_targ_bin = zeros(nSess,nBlocksSess,nBinsTarg);
    count_resp_bin = zeros(nSess,nBlocksSess,nBinsTarg);
    % also making sure that each individual position is used precisely
    % once per session. 
    count_targ_stim = zeros(nSess, length(orient_space_targ));
    count_resp_stim = zeros(nSess, length(orient_space_targ));
    
    for se = 1:nSess
        for rr = 1:nBlocksSess
            for tt = 1:nTrialsBlock
                % counting to make sure each thing happens once
                targ = trialInfo.TargPosition(se,rr,tt);
                [~,targbin] = find(targ==binListTarg);
                count_targ_bin(se,rr,targbin) = count_targ_bin(se,rr,targbin) + 1;
                [~,targind] = find(targ==orient_space_targ);
                count_targ_stim(se,targind) = count_targ_stim(se,targind) + 1;
                
                resp = trialInfo.ReportStartPosition(se,rr,tt);
                [~,respbin] = find(resp==binListTarg);
                count_resp_bin(se,rr,respbin) = count_resp_bin(se,rr,respbin) + 1;
                [~,respind] = find(resp==orient_space_targ);
                count_resp_stim(se,respind) = count_resp_stim(se,respind) + 1;
                
            end
        end
    end
    assert(all(count_targ_bin(:)==1));
    assert(all(count_resp_bin(:)==1));
    assert(all(count_targ_stim(:)==1));
    assert(all(count_resp_stim(:)==1));
    fprintf('All counterbalancing is good!\n')
    
    %% make plots to visualize things if we wish to
    if makePlots
        
        alltargpos = trialInfo.TargPosition;

        r = 1;

        figure;hold all;

        hist(alltargpos,360)
        title('target position distribution');
        xlim([0,360])

        set(gcf,'Color','w');
          
    end
    
    %% other info about timing

    t.ITIrange = [1,5];
    
    iti_list = linspace(t.ITIrange(1),t.ITIrange(2), nTrialsBlock);
    trialInfo.ITI = zeros(nSess,nBlocksSess,nTrialsBlock);
    for ss = 1:nSess
        for rr = 1:nBlocksSess
            iti_list_shuff = iti_list(randperm(nTrialsBlock));
            trialInfo.ITI(ss,rr,:) = iti_list_shuff;
        end
    end
    
    %% calculate the total length of each run (different for each b/c delay periods are different)
    % all this timing info will go straight into the task script, just to
    % make sure we get the right number of TRs
    
    t.PreTargetCueTime = 0.750; 
    t.TargetTime = 0.500;   
    t.DelayLength = 12;    
    t.MaxRespTime = 3;  % continuous report task time limit

    t.BeginFixation = 13;
    t.EndFixation = 8;

    % total run length
    t.trialLength = t.PreTargetCueTime + t.TargetTime+...
        t.DelayLength + t.MaxRespTime + mean(t.ITIrange);
    t.totalTime = t.BeginFixation + t.trialLength*nTrialsBlock + t.EndFixation;

    t.TR = 0.800;  
    t.totalTRsEachRun = ceil(t.totalTime/t.TR);

    %% save in structure to load
    
    trialInfo.t = t;
    
end

