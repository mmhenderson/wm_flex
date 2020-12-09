function trialInfo = OriSpin2_CounterBal_FMRI(subNum,makePlots)
% make a counterbalanced sequence of blocks over single session of 
% OriSpin2 (spatial WM) task

% trialInfo is a structure array [nBlocksTot x 1]

% each run has a unique task condition that is in effect for entire
% run.

% over all sessions counterbalanced for
    % response line position 
    % response mapping
    % distribution of memory stimuli (positions)

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

    datadir_local = [pwd filesep 'Data' filesep];
    fnsave = [datadir_local 'OriSpin2_S', subNum '_', trialInfo.dateStr '_MainTask.mat'];
    if exist(fnsave,'file')
%         error('cannot overwrite an existing file')
    end

    dbstop if error
    %% define some parameters
    nSess = 1;
    nBlocksTotal = 10;  % how many blocks to split trials into?
    nTasks = 2;
    nBinsTarg = 10; % for counter-balancing, how mamy bins do we divide the stim space into?
    nBinsBound = 10; % how many bins is boundary space divided into?
    nDegrees = 360;
    % how many degrees apart must the possible positions be? want this to
    % be small, below JND so subjects don't learn the grid
    axis_spacing = nDegrees/nBinsTarg/nBinsBound;
    
    % boundaries and orientations are staggered by a tiny bit, so that we
    % never have an ambiguous trial.
    orient_space_targ = axis_spacing:axis_spacing:nDegrees;
    orient_space_bound = axis_spacing/2:axis_spacing:nDegrees;
    
    binListTarg = reshape(orient_space_targ, ...
        length(orient_space_targ)/nBinsTarg, nBinsTarg);
    binListBound = reshape(orient_space_bound, ...
        length(orient_space_bound)/nBinsBound, nBinsBound);
    
    %% calculate more numbers
    nBlocksSess = nBlocksTotal/nSess;
    nTrialsTotal = nBinsTarg*nBinsBound*nTasks;
    nTrialsBlock = nTrialsTotal/nBlocksTotal;
    nTrialsEachTask = nTrialsTotal/nTasks;
    
    trialInfo.nSess = nSess;
    trialInfo.nTrialsBlock = nTrialsBlock;
    trialInfo.nBlocksSess = nBlocksSess;
    trialInfo.nTasks = nTasks;
    trialInfo.binListBound = binListBound;
    
    %% do all the counter-balancing here
    % initalize these things over all trials at once
    bound_position = zeros(nSess, nBlocksSess, nTrialsBlock);
    targ_position = zeros(nSess,nBlocksSess, nTrialsBlock);
    correct_resp = zeros(nSess,nBlocksSess, nTrialsBlock);

    % defining my order of tasks over all sessions/blocks
    % (random interleaved trials)
    taskinds = repelem([1:nTasks], nTrialsEachTask);
    taskinds_rand = taskinds(randperm(nTrialsTotal));
    
    % now looping over tasks
    for task = 1:nTasks
        
        resp_balanced = 0;
        maxiter = 10000;
        it=0;
        while ~resp_balanced && it<maxiter
            it=it+1;
            % first take my binned list of boundary positions - shuffle them
            % within bins (within columns). Result is [nBinsTarg x nBinsBound]
            bound_pos_this_task = binListBound;
            for bb = 1:nBinsBound
               bound_pos_this_task(:,bb) = binListBound(randperm(nBinsTarg),bb);
            end

            % now take my binned list of target positions, transpose it, and 
            % shuffle the rows. Result is [nBinsTarg x nBinsBound]
            targ_pos_this_task = binListTarg';
            for tt = 1:nBinsTarg
               targ_pos_this_task(tt,:) = binListTarg(randperm(nBinsBound),tt);
            end

            % shuffling over all trials in this task (over all sessions and
            % blocks and trials)
            trial_order_rand = randperm(nTrialsEachTask);
            bound_pos_this_task = reshape(bound_pos_this_task(trial_order_rand),size(bound_pos_this_task));
            targ_pos_this_task = reshape(targ_pos_this_task(trial_order_rand),size(targ_pos_this_task));

            % calculate what correct response should be
            correct_resp_this_task = zeros(size(targ_pos_this_task));
            over180 = bound_pos_this_task>180;
            correct_resp_this_task(over180 &...
                targ_pos_this_task<bound_pos_this_task &...
                targ_pos_this_task>mod(bound_pos_this_task+180,360)) = 2;
            correct_resp_this_task(over180 & ...
                (targ_pos_this_task>bound_pos_this_task |...
                targ_pos_this_task<mod(bound_pos_this_task+180,360))) = 1;
            correct_resp_this_task(~over180 &...
                targ_pos_this_task>bound_pos_this_task &...
                targ_pos_this_task<mod(bound_pos_this_task+180,360)) = 1;
            correct_resp_this_task(~over180 &...
                (targ_pos_this_task<bound_pos_this_task |...
                targ_pos_this_task>mod(bound_pos_this_task+180,360))) = 2;
            resp_balanced = sum(correct_resp_this_task(:)==1)==(nBinsBound*nBinsTarg/2)...
                & sum(correct_resp_this_task(:)==2)==(nBinsBound*nBinsTarg/2);
        end
        if it==maxiter
            error('failed to find a sequence that works!!')
        else
            fprintf('Task %d: responses balanced on iteration %d!\n',task,it)
        end
        %put them into the big array where they are supposed to be
        bound_position(taskinds_rand==task) = bound_pos_this_task;
        targ_position(taskinds_rand==task) = targ_pos_this_task;
        correct_resp(taskinds_rand==task) = correct_resp_this_task;
    end

    % these are the final variables we will reference during the experiment
    trialInfo.TargPosition = targ_position;
    trialInfo.BoundPosition = bound_position;
    trialInfo.Task = reshape(taskinds_rand,size(targ_position));
    trialInfo.CorrectResponse = correct_resp;
    
    %% figure out which trials lie in the extra-hard category 
    % give bonus $$ for getting these right!
    is_hard = zeros(size(targ_position));
    for se = 1:nSess
        for rr = 1:nBlocksSess
            for tt = 1:nTrialsBlock
                targ = trialInfo.TargPosition(se,rr,tt);
                bound = trialInfo.BoundPosition(se,rr,tt);
                [~,boundbin] = find(bound==binListBound);
                [~,targbin] = find(targ==binListTarg);
                is_hard(se,rr,tt) = boundbin==targbin;
            end
        end
    end
    assert(sum(is_hard(:))==nTasks*nBinsTarg);
    trialInfo.BonusTrial = is_hard;
    %% double check my counter-balancing

    count = zeros(nTasks, nBinsTarg, nBinsBound);
    
    for se = 1:nSess
        for rr = 1:nBlocksSess
            for tt = 1:nTrialsBlock
                % counting to make sure each combination happens once
                targ = trialInfo.TargPosition(se,rr,tt);
                bound = trialInfo.BoundPosition(se,rr,tt);
                [~,boundbin] = find(bound==binListBound);
                [~,targbin] = find(targ==binListTarg);
                task = trialInfo.Task(se,rr,tt);
                count(task,targbin,boundbin) = count(task,targbin,boundbin) + 1;
                % also check response expected
                if trialInfo.CorrectResponse(se,rr,tt)==1
                    if bound>180
                        assert(targ>bound || targ<mod(bound+180,360));
                    else
                        assert(targ>bound && targ<mod(bound+180,360));
                    end
                elseif trialInfo.CorrectResponse(se,rr,tt)==2
                    if bound>180
                        assert(targ<bound && targ>mod(bound+180,360));
                    else
                        assert(targ<bound || targ>mod(bound+180,360));
                    end
                end
            end
        end
    end
    assert(all(count(:)==1));
    fprintf('All counterbalancing is good!\n')
    
    %% make plots to visualize things if we wish to
    if makePlots
        for task=1:nTasks
            alltargpos = trialInfo.TargPosition(trialInfo.Task==task);
            allboundpos = trialInfo.BoundPosition(trialInfo.Task==task);
            r = 1;
            
            figure;hold all;
            subplot(1,2,1);hold all;
            hist(alltargpos,360)
            title('target position distribution');
            xlim([0,360])
             
            subplot(1,2,2);hold all;
            hist(allboundpos,360)
            title('boundary position distribution');
            xlim([0,360])
            suptitle(sprintf('Task %d',task));
            set(gcf,'Color','w');
        end       
    end
    
 
    %% another plot (scatter)
    if makePlots
        axlims = [-5, 365];
        for task=1:nTasks
            figure();hold all;
            alltargpos = trialInfo.TargPosition(trialInfo.Task==task);
            allboundpos = trialInfo.BoundPosition(trialInfo.Task==task);            
            plot(alltargpos,allboundpos,'.')
            for tt = 1:nBinsTarg
                plot([binListTarg(1,tt)-axis_spacing/2,binListTarg(1,tt)-axis_spacing/2],...
                   axlims,'-','Color','k');
            end
            plot([binListTarg(end,end)+axis_spacing/2,binListTarg(end,end)+axis_spacing/2],...
                axlims,'-','Color','k');
            for bb = 1:nBinsBound
                plot(axlims,[binListBound(1,bb)-axis_spacing/2,binListBound(1,bb)-axis_spacing/2],...
                    '-','Color','k');
            end
            plot(axlims,[binListBound(end,end)+axis_spacing/2,binListBound(end,end)+axis_spacing/2],...
                '-','Color','k');

            axis square
            xlim(axlims)
            ylim(axlims)
            xlabel('Spatial memory positions (degrees)');
            ylabel('Response boundary positions (degrees)');
            title(sprintf('Task %d\nGridded, counterbalanced stimulus positions',task))
            set(gcf,'Color','w');
        end
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
    t.PreCueDelayLength = 1;
    t.CueLength = 2;
    t.BoundPreviewLength = 1;
    t.PostCueDelayLength = 12;
    t.BoundFinalLength = 2;
    t.MaxRespTime = 3;  % end of this happens during ITI

    assert(min(t.ITIrange)>=(t.MaxRespTime-t.BoundFinalLength))
    
    t.BeginFixation = 13;
    t.EndFixation = 8;

    % total run length
    t.trialLength = t.PreTargetCueTime + t.TargetTime+t.PreCueDelayLength+t.CueLength+...
        t.BoundPreviewLength+t.PostCueDelayLength+t.BoundFinalLength + mean(t.ITIrange);
    t.totalTime = t.BeginFixation + t.trialLength*nTrialsBlock + t.EndFixation;

    t.TR = 0.800;  
    t.totalTRsEachRun = ceil(t.totalTime/t.TR);

    %% save in structure to load
    
    trialInfo.t = t;

    save(fnsave,'trialInfo');
    
end

