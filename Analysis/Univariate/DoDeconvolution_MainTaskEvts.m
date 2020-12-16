%% Deconvolution on main task data. 
% estimate BOLD resp amplitude for each type of event in each condition. 

clear
close all;
sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
nROIs_all = numel(ROI_names);

plotVisMotor = 1;
plotVisMotorSS = 1;

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));

nROIs = length(plot_order_all);

col = plasma(5);
col = col(2:2:end-1,:);
lw=1;

condLabStrs = {'Predictable','Random'};

% length of each run type
nTRs = 583-16;
trDur = 0.80;

% deconvolution param
nTRs_out = 30; % # of TRs post-stimulus to model in the decon analysis
tax = trDur*(0:nTRs_out-1); 

nTrialsPerRun = 20;
nConds = 2;

allsub_HRFs_mean = zeros(nSubj, length(plot_order_all), nConds, nTRs_out);
allsub_HRFs_sem = zeros(nSubj, length(plot_order_all), nConds, nTRs_out);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.62,0.655,0.69];

col_conds = [125, 93, 175; 15, 127, 98]./255;

diff_col=[0.5, 0.5, 0.5];

ylims=[-0.7, 0.7];
fs=20;

%%
for ss=1:length(sublist)

    
    substr = sprintf('S%02d',sublist(ss));
   
    fn2load = fullfile(exp_path,'Samples',sprintf('SampleFile_%s.mat',substr));
    load(fn2load, 'samplesMain','main','ROIs','all_vox_concat');
    
    %% load the timing file (made in GetEventTiming.m)
    
    fn = fullfile(exp_path,'Samples',sprintf('TimingFile_%s.mat',substr));
    if ~exist(fn, 'file')
        error('need to make timing file first, run GetEventTiming.m')
    end
    fprintf('Loading event timing file\n')
    load(fn)

    %%
    for vi = 1:nROIs % for all visual areas I want to look at
        vv = plot_order_all(vi);
        if strcmp(ROI_names{vv},'PMc')
            name2use='Premotor';
        else 
            name2use=ROI_names{vv};
        end
        %% pull out the data from each ROI
        % want both hemispheres
        [rowind1,colind1] = find(strcmp(reshape({ROIs.name},2,nROIs_all),sprintf('lh_%s',name2use)));
        [rowind2,colind2] = find(strcmp(reshape({ROIs.name},2,nROIs_all),sprintf('rh_%s',name2use)));
        col_inds = [colind1,colind2]; % column is the region
        row_inds = [rowind1,rowind2];   % row is the hemisphere
        
        if vv<12
            % make sure each early visual area is defined for both
            % hemispheres.
            assert(numel(col_inds)==2)
        end
        
        mainDat=[]; 
        for ii=1:length(col_inds)
            name = ROIs(row_inds(ii),col_inds(ii)).name;
            if ~isempty(ROIs(row_inds(ii),col_inds(ii)).voxel_inds)
                % jj gives indices into the all_vox_concat array
                [~,jj]=intersect(all_vox_concat, ROIs(row_inds(ii),col_inds(ii)).voxel_inds);
                mainDat = [mainDat, samplesMain(:,jj)];
            end
        end
        nVox = size(mainDat,2);

        if nVox==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end

        fprintf('processing area %s, %d voxels\n', ROI_names{vv}, nVox);
        
        %% now zscore the data from each run to normalize...
        
        nRuns = size(mainDat,1)/nTRs; % hopefully 5 runs per session
        if mod(nRuns,1)~=0
            error('something bad happened here with mainDat run length')
        end
        for ii=1:nRuns
            mainDat(ii*nTRs-nTRs+1:ii*nTRs,:) = zscore(mainDat(ii*nTRs-nTRs+1:ii*nTRs, :),1);
        end
        
        assert(numel(unique(main.RunLabels))==nRuns);
        %% label the data
        % event labels are as follows:
        % pre-targ cue, targ, delay1, cue, bound-preview, delay2,
        % bound-actual, start ITI
        % Event_type = [Event_type, [0.2, 1, 0, 0.3, 2, 0, 3, 0]];
        event_labels_reshaped = reshape(main.EventLabels,nTRs,length(main.EventLabels)/nTRs);

        % now find the actual onset of each trial - switch from 0.2 to 1
        % (or 0 to 1)
        trial_onset_bool = event_labels_reshaped==1;
        trial_onset_bool = trial_onset_bool(:);
        trial_onset_num = find(trial_onset_bool);

        nTrials = nRuns*nTrialsPerRun;
        assert(numel(trial_onset_num)==nTrials);
        
        %% get conditions

        wmConds = zeros(size(main.EventLabels));
        wmConds(trial_onset_num) = main.CondLabels(trial_onset_num); 
        
        unconds = unique(wmConds(wmConds~=0));
        for uu=1:numel(unconds)
            wmConds(wmConds==unconds(uu)) = uu;
        end

        condLabStrs = {'Targ Onset: Predictable','Targ Onset: Random'};
        %% now do deconvolution to estimate the event related HRF associated with each event type.
        % written in a very general way that can be applied to any data set -
        % i.e. will be explicit about computing #of conditions, and will also
        % include a constant term for each run even though we've zero-meaned 
        % the data (so this model will be appropriate for analyzing raw data 
        % that has large baseline shifts between runs as well).

        nConds = numel(unique(wmConds))-1;      % minus 1 because we're not modelling 0's, just donuts and small circles.
        nMainRuns = numel(unique(main.RunLabels));   % can do this size(oriLocDat,1)/oriLocTRs, but we'll use this method so that we can cross-check each approach (see error checking inside of the doDecon func)

        maintaskHRFs = doDecon_ForRR(mainDat, wmConds, nConds, nMainRuns, nTRs, nTRs_out);
        for cc=1:nConds
            allsub_HRFs_mean(ss,vi,cc,:) = mean(maintaskHRFs(:,cc,:),3);
            allsub_HRFs_sem(ss,vi,cc,:) = std(maintaskHRFs(:,cc,:),[],3)./sqrt(nVox);
        end
        
       
    end

end

%% compare between conditions
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 122434;
rng(rndseed,'twister')
nPermIter=1000;
vals = allsub_HRFs_mean;
real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);

% p_diff_sr=nan(nROIs,nTRs_out);
for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,tr));
        % what is the sign-rank statistic for the real data?
        real_sr_stat(vv,tr) = signrank_MMH(realvals(:,1),realvals(:,2));
        % determine before the parfor loop which conditions get randomly
        % swapped on each iteration (otherwise not deterministic)
        inds2swap = double(randn(nSubj,nPermIter)>0);
        inds2swap(inds2swap==0) = -1;
        parfor ii=1:nPermIter          
            % randomly permute the condition labels within subject
            randvals=realvals;
            randvals(inds2swap(:,ii)==-1,:) = randvals(inds2swap(:,ii)==-1,[2,1]);    
            % what is the sign-rank statistic for this randomly permuted data?
            rand_sr_stat(vv,tr,ii) = signrank_MMH(randvals(:,1),randvals(:,2));
        end
    end
end

% % compute a two-tailed p-value comparing the real stat to the random
% % distribution. Note that the <= and >= are inclusive, because any
% % iterations where real==null should count toward the null hypothesis. 
p_cond_diff = 2*min(cat(3,mean(repmat(real_sr_stat,1,1,nPermIter)>=rand_sr_stat,3), ...
    mean(repmat(real_sr_stat,1,1,nPermIter)<=rand_sr_stat,3)),[],3);

array2table(squeeze(p_cond_diff(vismotor_inds,1:15)),...
    'RowNames',vismotor_names,'VariableNames',strseq('CondDiff_TR',1:15))
array2table(squeeze(p_cond_diff(vismotor_inds,16:30)),...
    'RowNames',vismotor_names,'VariableNames',strseq('CondDiff_TR',16:30))

%%
if plotVisMotor
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));
    figure;hold all;
      
    for vi = 1:numel(vismotor_inds)
        subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
        lh=[];lx=0;ll=[];
        for cc=1:nConds
          
            % take contra - ipsi difference
            vals = allsub_HRFs_mean(:,vismotor_inds(vi),cc,:);

            if nSubj==1
                meanvals =squeeze(vals);
                semvals = squeeze(allsub_HRFs_sem(:,vismotor_inds(vi),cc,:));
            else
                meanvals = squeeze(nanmean(vals,1));
                semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
            end

            lh=[lh,plot(tax,meanvals,'-','Color',col_conds(cc,:),'LineWidth',lw)];
            bandedError_MMH(tax, meanvals',semvals', col_conds(cc,:), 0.5);
%                 lh=[lh, errorbar(tax, meanvals, semvals,'Color',col_conds(cc,:,bb),'LineWidth',lw)];
            lx=lx+1;
            ll{lx} = sprintf('%s',condLabStrs{cc});
            
        end
        
        for aa=1:numel(alpha_vals)
            inds2plot=p_cond_diff(vismotor_inds(vi),:)<alpha_vals(aa);
            plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa));
            lh=[lh, plot(tax(1)-5, sig_heights(nConds+1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
            ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
        end
                
            
        set(gca, 'FontSize', 12, 'XLim',[0 max(tax)],'YLim',ylims)
        h=plot(get(gca,'XLim'),[0,0],'--','Color',[0.7, 0.7, 0.7]);
        uistack(h,'bottom')
        for ee = 1:length(evts2plot)
            h=plot([evts2plot(ee),evts2plot(ee)],ylims,'-','Color',[0.8, 0.8, 0.8]);
            uistack(h,'bottom')
        end
        h=fill([evts2plot(1),evts2plot(2), evts2plot(2),evts2plot(1)], repelem(ylims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        h=fill([evts2plot(3),evts2plot(4), evts2plot(4),evts2plot(3)], repelem(ylims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        if vi==1
            xlabel('Time(s)');
            ylabel('BOLD resp');
        end

        if vi==numel(vismotor_inds)
            legend(lh, ll);
        end

        
        title(sprintf('%s',vismotor_names{vi}));
        
        
    end
    suptitle('Deconvolved BOLD (HRF) each condition')
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200]);
    
    
end
%%
subcolors = viridis(nSubj+1);
sslims =  [-0.5, 0.9];
if plotVisMotorSS
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));

    for cc=1:nConds
        figure;hold all;
        for vi = 1:numel(vismotor_inds)
            subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
            ll=[];lx=0;lh=[];
            for ss=1:nSubj
              
                vals = squeeze(allsub_HRFs_mean(ss,vismotor_inds(vi),cc,:));
                
                lh=[lh,plot(tax,vals,'-','Color',subcolors(ss,:),'LineWidth',lw)];
                ll{ss}=sprintf('S%02d',sublist(ss));
                
            end
                
            set(gca, 'FontSize', 12, 'XLim',[0 max(tax)],'YLim',sslims)
            h=plot(get(gca,'XLim'),[0,0],'--','Color',[0.7, 0.7, 0.7]);
            uistack(h,'bottom')
            for ee = 1:length(evts2plot)
                h=plot([evts2plot(ee),evts2plot(ee)],sslims,'-','Color',[0.8, 0.8, 0.8]);
                uistack(h,'bottom')
            end
            h=fill([evts2plot(1),evts2plot(2), evts2plot(2),evts2plot(1)], repelem(sslims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
            uistack(h,'bottom')
            h=fill([evts2plot(3),evts2plot(4), evts2plot(4),evts2plot(3)], repelem(sslims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
            uistack(h,'bottom')
            if vi==1
                xlabel('Time(s)');
                ylabel('BOLD resp');
            end

            if vi==numel(vismotor_inds)
                legend(lh,ll);
            end

            if contains(vismotor_names{vi}, ' ')
                % break it into two strings
                spaceind = find(vismotor_names{vi}==' ');
                title(sprintf('%s\n%s', vismotor_names{vi}(1:spaceind-1), vismotor_names{vi}(spaceind+1:end)));
            else
                title(sprintf('%s',vismotor_names{vi}));
            end
        end
        suptitle(sprintf('%s - Deconvolved BOLD (HRF)',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1800,1200]);
    end
end