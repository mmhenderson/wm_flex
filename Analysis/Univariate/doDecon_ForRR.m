function hrfs = doDecon_ForRR(dat, condVec, nConds, nRuns, nTRs, nTRsToModel) 
%
% This function does deconvolution of functional MRI data to estimate the
% Hemodynamic Response Function (HRF) associated with each event type.
%
% Input:
% dat --> data matrix (time x voxel). i.e. one voxel's timeseries in each 
% column (so time goes down the rows!)
% condVec --> vector of time x 1 (i.e. size(dat,1),1) filled with 0's, with 
% non-zero entries in cells representing the onset TR of each event and the
% non-zero value will indicate the event type (should be integers from
% 1:nConds for this function to work)
% nConds --> number of conditions to model (== to unique(condVec)-1) 
% nRuns --> number of runs, each having nTRs (can figure this out in this
% function, but may as well be explicit in passing)
% nTRs --> number of TRs in each run
% nTRsToModel --> number of TRs to model after event onset
%
% written js 03082017 for RR's WM data. RR doing super bug testing on
% 03132017 to make sure it can be used without care.

% do a few quick checks
if size(dat,1)/nTRs ~= nRuns
    fprintf('\nWrong nTRs/nRuns/dat passed into doDecon_ForRR.m\nReturning...\n')
    return;
end

if sum(unique(condVec)~=0) ~= nConds % all unique conditions in condVec that are not 0 
    fprintf('\nbad condVec (or wrong nCond) passed into doDecon_ForRR.m\nReturning...\n')
    return;
end

if numel(condVec) ~= size(dat,1)
    fprintf('\ndimensions of condVec and dat don''t match\nReturning...\n')
    return;
end

% # of voxels
nVox = size(dat,2);


%% do deconvolution to estimate the event related HRF associated with each event type.
% going to try to write in a very general way that can be applied to any
% data set - i.e. will be explicit about computing #of condtions, and will
% include a constant term for each run even though we've zero-meaned the
% data (so this model would be appropriate for analyzing raw data that has
% large baseline shifts between runs as well). 


%% make a design matrix
% allocate a DM (Design Matrix) "X" with enough columns for the constant 
% terms. By convention the first columns will correspond to the shifted 
% timeseries for each condition (i.e. the regressors to estimate the HRF 
% for each condition) and then the last nRuns columns will contain the 
% constant term for each run). Note that the constant terms are separate 
% for each run to support the use of raw data where there is a mean offset
% between runs --> even though this mean offset is removed by z-scoring or 
% PSC (Percent Signal Change) conversion, it does not hurt to have these 
% constants in the model still (they just turn out to be really really 
% really close to zero... and if they're not, then you know something went 
% wrong!). 
X = zeros(size(dat,1), nConds*nTRsToModel+nRuns); 

% then we do a triple loop to fill up the DM - the outer loops over runs, 
% then the next loops over conditions, and the inner loops over time-points 
% in our modelled HRF
for r = 1:nRuns
    
    % define an index into the timepoints that correspond to the current run
    runInd = r*nTRs-nTRs+1:r*nTRs;
    
    % first put in the constant term for the model
    X(runInd, nConds*nTRsToModel+r) = ones(nTRs,1);
    
    % then loop over conditions and put in a shifted version of each
    % condition's timeseries in the appropriate column
    for c = 1:nConds
        
        % get a column vector of 0's and 1's marking the onset time of the
        % all events in the current condition
        evtVec = zeros(nTRs,1); % make a blank tmp column vector to store trial sequence from this run...
        evtVec(condVec(runInd)==c) = 1; % insert a 1 to mark the onset of all instances of the current condition
                
        % then loop over the number of TRs in the HRF that you want to model
        for t = 1:nTRsToModel
            
            % insert the temporally shifted event timeseries from this run.
            % IMPORTANT: make sure that the shifted timeseries DO NOT 
            % extend beyond the end of the run OR wrap around into itself 
            % to the begining of the run, or else you'll model the last 
            % events as causing variance in the BOLD signal at the begining 
            % of the run. 
            X(runInd,t+(c*nTRsToModel-nTRsToModel)) = [zeros(t-1,1); evtVec(1:end-(t-1))]; % note that we're just cutting 1 entry off the end of the vector each time (and padding zeros in front to keep vector length constant). 
        end
    end
end
% DM build done...


% check to make sure that nothing super bad happened...always look at it
% first: imagesc(X) will show the whole thing, including per-run regressors
% second: imagesc(X(1:175,:)) will show a zoomed in version in time. Should
% look orderly as you'll see... especially check out the transition zone 
% between two runs to make sure that the regressors stop at the of a run 
% and don't continue on past the end of a run. 
% third: imagesc(X(nTRs*2-nTRsToModel*2:nTRs*2+nTRsToModel*2,:)) will allow
% you to check that the constant term shifts over 1 column at the end of 
% each run (although you can also look at imagesc(X) and see that there are 
% nRuns vertical bars at the right side ofthe figure). 
% last: Take a last take a look at the DM for a few runs like this:
% imagesc(X(nTRs:2*nTRs,:)) or imagesc(X(9*nTRs:10*nTRs,:))
% After plotting, do at least 1 check to make sure that its not rank
% deficient (i.e. mae sure rank == number of columns or regressors).
if rank(X) ~= size(X,2)
    display('Rank Deficient design matrix for the WM runs')
    return;
end


%% then solve the GLM 
% and reshape to get the estimated event-related HRF for each condition 
% (and note that this will compute the model for every voxel as well). 
b = X\dat; % beta weight estimated for each timepoint and condition
% to illustrate what this does: solving this per voxel would look like
% voxel_b = X\dat(:,voxel) - where voxel_b will give the estimated beta
% weights for that single voxel based off the entire DM "X" and the data
% for only that one voxel (i.e. dat(:,voxel)). So that X multiplied with
% the voxel_b beta weights will give you the activity over time in that
% voxel dat(:,voxel) --> i.e. b*X = dat


% we just care about the first nConds * nTRsToModel beta weights cause the
% last nRuns weights are just the constant terms for each run. 
hrfs = NaN(nTRsToModel, nConds, nVox); % to store event related data from each condition and voxel
for c = 1:nConds
    curCondInd = c*nTRsToModel-nTRsToModel+1:c*nTRsToModel; % indices for the current condition
    hrfs(:,c,:) = b(curCondInd,:);
end
% check that ~sum(isnan(hrfs(:)))

% check the range of my offset parameters
% figure; hist(mean(b(curCondInd(end)+1:end,:)));
constant_range = [min(min(b(curCondInd(end)+1:end,:))) max(max(b(curCondInd(end)+1:end,:)))];
hfrs_range = [min(min(min(hrfs))) max(max(max(hrfs)))];




