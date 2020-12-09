% extract signal in each voxel averaging over several TRs following trial
% events, save values for input to further analysis

clear
close all

sublist = [7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

%% params here
trDur = .8; % actual duration of the TR, for plotting, in seconds...
nTRs = 452-16;

ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','Premotor'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1-S1 all'};


nROIs = length(ROI_names);

% which TRs am i averaging over? From the target onset time.
avgTRs_stim = [3,7];
% also averaging over data after the response probe onset - get a response
% related HRF in motor areas.
avgTRs_resp = [4,7];

nTRs2concat = 25;

nTrialsPerRun=20;
    
%% load data    
for ss = sublist

    locSig = struct();

    substr = sprintf('S%02d',ss);
   
    fn2load = fullfile(exp_path,'Samples',sprintf('SampleFile_%s.mat',substr));
    load(fn2load, 'samplesDWM','ROIs','all_vox_concat');
    
    %% load the timing file (made in GetEventTiming.m)
    
    fn = fullfile(exp_path,'Samples',sprintf('TimingFile_%s.mat',substr));
    if ~exist(fn, 'file')
        error('need to make timing file first, run GetEventTiming.m')
    end
    fprintf('Loading event timing file\n')
    load(fn,'dwm')
    
    fn2save = fullfile(exp_path,'Samples',sprintf('DWMLocSignalByTrial_%s.mat',substr));    
   
    %%
    for vv = 1:nROIs % for all visual areas I want to look at
        %% pull out the data from each ROI
        % want both hemispheres
        [rowind1,colind1] = find(strcmp(reshape({ROIs.name},2,nROIs-1),sprintf('lh_%s',ROI_names{vv})));
        [rowind2,colind2] = find(strcmp(reshape({ROIs.name},2,nROIs-1),sprintf('rh_%s',ROI_names{vv})));
        col_inds = [colind1,colind2]; % column is the region
        row_inds = [rowind1,rowind2];   % row is the hemisphere
        
        if strcmp(ROI_names{vv},'M1-S1 all')
            inds=[12:14];
            row_inds=[];col_inds=[];
            for ii = inds

                [r1,c1] = find(reshape(contains({ROIs.name}, sprintf('lh_%s',ROI_names{ii})),2,[]));
                [r2,c2] = find(reshape(contains({ROIs.name}, sprintf('rh_%s',ROI_names{ii})),2,[]));
                if isempty(r1) & isempty(r2)
                    fprintf('no voxels in %s for %s\n',ROI_names{ii}, substr);
                end
                row_inds=[row_inds; r1;r2];
                col_inds=[col_inds; c1;c2];
            end
        end
        
        if vv<12
            % make sure each early visual area is defined for both
            % hemispheres.
            assert(numel(col_inds)==2)
        end
        
        locDat=[]; 
%         eccVals = [];
        for ii=1:length(col_inds)
            name = ROIs(row_inds(ii),col_inds(ii)).name;
            if ~isempty(ROIs(row_inds(ii),col_inds(ii)).voxel_inds)
                % jj gives indices into the all_vox_concat array
                [~,jj]=intersect(all_vox_concat, ROIs(row_inds(ii),col_inds(ii)).voxel_inds);
                locDat = [locDat, samplesDWM(:,jj)];
            end
        end
        nVox = size(locDat,2);

        if nVox==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end

        fprintf('processing area %s, %d voxels\n', ROI_names{vv}, nVox);
        
        %% now zscore the data from each run to normalize...

        nRuns = size(locDat,1)/nTRs; % hopefully 5 runs per session
        if mod(nRuns,1)~=0
            error('something bad happened here with DWM run length')
        end
        for ii=1:nRuns
            locDat(ii*nTRs-nTRs+1:ii*nTRs,:) = zscore(locDat(ii*nTRs-nTRs+1:ii*nTRs, :));
        end
        
        assert(numel(unique(dwm.RunLabels))==nRuns);
        %% label the data
        % event labels are [1 for target, 0 for delay, 3
        % for response, 0 for ITI].
        trial_labels_reshaped = reshape(dwm.TrialLabels,nTRs,length(dwm.TrialLabels)/nTRs);
        
        % now find the actual onset of each trial - switch from 0.2 to 1
        % (or 0 to 1)
        trial_onset_bool = [ones(1,size(trial_labels_reshaped,2));diff(trial_labels_reshaped)==1];
        trial_onset_bool = trial_onset_bool(:);
        trial_onset_num = find(trial_onset_bool);
        
        event_labels_reshaped = reshape(dwm.EventLabels,nTRs,length(dwm.EventLabels)/nTRs);

        resp_onset_bool = [zeros(1, size(event_labels_reshaped,2)); diff(event_labels_reshaped)==3];
        resp_onset_bool = resp_onset_bool(:);
        resp_onset_num = find(resp_onset_bool);
        
        nTrials = nRuns*nTrialsPerRun;
        assert(numel(trial_onset_num)==nTrials);
        assert(numel(resp_onset_num)==nTrials);
        %% save out a bunch of descriptors for the trials
        
        locSig(vv).runLabs = dwm.RunLabels(trial_onset_num);
        locSig(vv).trialLabs = dwm.TrialLabels(trial_onset_num);
        locSig(vv).ExpDigit = dwm.ExpDigit(trial_onset_num);
        locSig(vv).ActDigit = dwm.ActDigit(trial_onset_num);
        locSig(vv).RespEarly = dwm.RespEarly(trial_onset_num);
        locSig(vv).RT = dwm.RT(trial_onset_num);
        
        %% avg the data across each trial
        % single value for each trial, averaged over multiple TRs
        dat_avg_targ = nan(nTrials, nVox); 
        dat_avg_resp = nan(nTrials, nVox); 
       
        % TR-by-TR data        
        dat_by_TR = nan(nTrials, nTRs2concat, nVox);
        
        
        triCnt = 0; % counter across "good" trials where the entire desired avg window is available.
    
        for rr=1:nRuns

            runInd = rr*nTRs-nTRs+1:rr*nTRs;
            assert(all(find(dwm.RunLabels==rr)==runInd'))
            curDat = locDat(runInd,:);      % data from current run

            % get numeric indices for each event in this run
            these_targ_onsets = find(trial_onset_bool(runInd));       
            these_resp_onsets = find(resp_onset_bool(runInd));
            
            assert(numel(these_targ_onsets)==numel(these_resp_onsets))
            assert(numel(these_targ_onsets)==nTrialsPerRun);
            
            for tt=1:numel(these_targ_onsets)
                % check to make sure we don't go past the end of the run.
                if these_targ_onsets(tt)+avgTRs_stim(2)<=nTRs && these_resp_onsets(tt)+avgTRs_resp(2)<=nTRs
               
                    triCnt = triCnt + 1;  % increment counter over good trials

                    % sum the data over desired interval.
                    dat_avg_targ(triCnt,:) = mean(curDat(these_targ_onsets(tt)+avgTRs_stim(1):these_targ_onsets(tt)+avgTRs_stim(2),:));
                    dat_avg_resp(triCnt,:) = mean(curDat(these_resp_onsets(tt)+avgTRs_resp(1):these_resp_onsets(tt)+avgTRs_resp(2),:));
                    
                end
                
                % also collecting data at each timept
                for tr = 1:nTRs2concat
                    if these_targ_onsets(tt)+tr-1<=nTRs 
                        dat_by_TR(triCnt,tr,:) = curDat(these_targ_onsets(tt)+tr-1,:);
                    else
                        error('TR %d is past the end of your run!!', tr)
                    end
                end
                
            end
        end

        assert(triCnt==nTrials)
        assert(~sum(isnan(dat_avg_targ(:))) && ~sum(isnan(dat_by_TR(:))))
        assert(~sum(isnan(dat_avg_resp(:))))

        locSig(vv).dat_avg_targ = dat_avg_targ;
        locSig(vv).dat_avg_resp = dat_avg_resp;

        % this matrix is [nTrials x nTRs x nVox]
        locSig(vv).dat_by_TR = dat_by_TR;
    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'locSig','ROI_names');

end