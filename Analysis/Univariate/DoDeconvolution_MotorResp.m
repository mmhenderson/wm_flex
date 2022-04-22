%% Plot signal related to contralateral/ipsilateral finger presses
% Get the average univariate signal in each hemisphere ROI on trials
% where either the left or right finger was to be pressed at end of delay. 
% See when/where univariate motor signals emerge.

% Used to generate Figure 1 - figure supplement 2.
%%
clear
close all;
sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
figpath = fullfile(exp_path,'figs');
addpath(fullfile(exp_path,'Analysis','stats_code'));
addpath(fullfile(exp_path,'Analysis','plotting_utils'));

% names of the ROIs - note the motor areas have been changed to nicer names
% here for plotting.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
nROIs_all = numel(ROI_names);

hemi_names={'lh','rh'};


plotVisMotor = 1;
plotVisMotorSS=1;
plotVisMotor_Diff = 1;
plot_order1 = [1:5,10,11,6:9,12:14];
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
nROIs = length(plot_order_all);

lw=1;

condLabStrs = {'Predictable','Random'};
respLabStrs = {'Contra finger','Ipsi finger'};

% length of each run type
nTRs = 583-16;
trDur = 0.80;

% deconvolution param
nTRs_out = 30; % # of TRs post-stimulus to model in the decon analysis
% define timing axis
tax = trDur*(0:nTRs_out-1); 

nTrialsPerRun = 20;
nConds = 2;
nResp=2;

% the nResp dimension goes [contra, ipsi, contra - ipsi]
allsub_HRFs_mean = nan(nSubj, length(plot_order_all), nConds, nResp+1, nTRs_out);
allsub_HRFs_sem = nan(nSubj, length(plot_order_all), nConds, nResp+1, nTRs_out);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.42,0.455,0.49];
diff_col=[0.5, 0.5, 0.5];

ylims = [-0.3, 0.7];
difflims = [-0.5 0.5];

col_resp = plasma(3);
col_resp = col_resp(1:2,:);

col_conds = [125, 93, 175; 15, 127, 98]./255;


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
        fprintf('processing %s area %s\n', substr,ROI_names{vv});
        if strcmp(ROI_names{vv},'PMc')
            name2use='Premotor';
        else 
            name2use=ROI_names{vv};
        end
        %% pull out the data from each ROI
        % want both hemispheres
        [rowind1,colind1] = find(strcmp(reshape({ROIs.name},2,nROIs_all),sprintf('lh_%s',name2use)));
        [rowind2,colind2] = find(strcmp(reshape({ROIs.name},2,nROIs_all),sprintf('rh_%s',name2use)));
       
        % jj gives indices into the all_vox_concat array
        [~,jj]=intersect(all_vox_concat, ROIs(rowind1,colind1).voxel_inds);
        mainDat_lh = samplesMain(:,jj);
        % jj gives indices into the all_vox_concat array
        [~,jj]=intersect(all_vox_concat, ROIs(rowind2,colind2).voxel_inds);
        mainDat_rh = samplesMain(:,jj);
      
        assert(size(mainDat_lh,2)>0 && size(mainDat_rh,2)>0)
        
        fprintf('processing area %s\n', ROI_names{vv});
        
        %% now zscore the data from each run to normalize...
        
        nRuns = size(mainDat_lh,1)/nTRs; % hopefully 5 runs per session
        if mod(nRuns,1)~=0
            error('something bad happened here with mainDat run length')
        end
        for ii=1:nRuns
            mainDat_lh(ii*nTRs-nTRs+1:ii*nTRs,:) = zscore(mainDat_lh(ii*nTRs-nTRs+1:ii*nTRs, :),1);
            mainDat_rh(ii*nTRs-nTRs+1:ii*nTRs,:) = zscore(mainDat_rh(ii*nTRs-nTRs+1:ii*nTRs, :),1);
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

        condLabs = zeros(size(main.EventLabels));
        condLabs(trial_onset_num) = main.CondLabels(trial_onset_num); 
        
        respLabs = zeros(size(main.CorrectResp));
        respLabs(trial_onset_num) = main.CorrectResp(trial_onset_num);
        
        unconds = unique(condLabs(condLabs~=0));
        unresp = unique(respLabs(respLabs~=0));
        
        % now making a big list w 4 conditions, for the two task conditions
        % and the left finger/right finger trials.
        hrfConds = zeros(size(respLabs));
        uu=0;
        cl=zeros(numel(unresp)*numel(unconds),1);
        rl=zeros(numel(unresp)*numel(unconds),1);
        for cc=1:numel(unconds)
            for rr=1:numel(unresp)
                uu=uu+1;
                cl(uu) = cc;
                rl(uu) = rr;
                hrfConds(condLabs==unconds(cc) & respLabs==unresp(rr)) = uu;
            end
        end

        %% now do deconvolution to estimate the event related HRF associated with each event type.
        % written in a very general way that can be applied to any data set -
        % i.e. will be explicit about computing #of conditions, and will also
        % include a constant term for each run even though we've zero-meaned 
        % the data (so this model will be appropriate for analyzing raw data 
        % that has large baseline shifts between runs as well).

        nHRFConds = numel(unique(hrfConds))-1;      % minus 1 because we're not modelling 0's, just donuts and small circles.
        nMainRuns = numel(unique(main.RunLabels));   % can do this size(oriLocDat,1)/oriLocTRs, but we'll use this method so that we can cross-check each approach (see error checking inside of the doDecon func)

        maintaskHRFs_lh = doDecon_ForRR(mainDat_lh, hrfConds, nHRFConds, nMainRuns, nTRs, nTRs_out);
        maintaskHRFs_rh = doDecon_ForRR(mainDat_rh, hrfConds, nHRFConds, nMainRuns, nTRs, nTRs_out);
        
        for cc=1:nConds
            for rr=1:nResp

                  % now switching from l/r finger labels to
                  % contralateral/.ipsilateral - so need to know which
                  % hemisphere to take which trials from.
                  
                  % rr=1 for contra
                  % want data from lh in brain, when right finger (respLabs==2) was used)
                  % want data from rh in brain, when left finger (respLabs==1) was used)
                  % rr=2 for ipsi
                  % want data from lh in brain, when left finger (respLabs==1) was used)
                  % want data from rh in brain, when right finger (respLabs==2) was used)
                  
                  uu_lh=find(cl==cc & rl==(3-rr));
                  uu_rh=find(cl==cc & rl==rr);
                  dat1 = squeeze(maintaskHRFs_lh(:,uu_lh,:));
                  dat2 = squeeze(maintaskHRFs_rh(:,uu_rh,:));
                  
                  hrf_concat = cat(2,dat1,dat2);
                  nVox = size(hrf_concat, 2);
                  
                  if rr==1
                      contra_hrf_concat = hrf_concat;
                  elseif rr==2
                      ipsi_hrf_concat = hrf_concat;
                  end
                  
                  allsub_HRFs_mean(ss,vi,cc,rr,:) = mean(hrf_concat,2);
                  allsub_HRFs_sem(ss,vi,cc,rr,:) = std(hrf_concat,[],2)./sqrt(nVox);
                  
            end
              
             % now take the contra-ipsi difference
            diffdat = contra_hrf_concat-ipsi_hrf_concat;

            % avg and sem over voxels
            allsub_HRFs_mean(ss,vi,cc,3,:) = mean(diffdat,2);
            allsub_HRFs_sem(ss,vi,cc,3,:) = std(diffdat,[],2)./sqrt(nVox);
        
        end
       
    end
end

%%
assert(~any(isnan(allsub_HRFs_mean(:))))
assert(~any(isnan(allsub_HRFs_sem(:))))

%% compare contra vs ipsi within each condition - wilcoxon signed rank tests
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 132434;
rng(rndseed,'twister')
nPermIter=1000;
real_sr_stat = nan(nROIs,nConds, nTRs_out);
rand_sr_stat = nan(nROIs, nConds, nTRs_out, nPermIter);
vals = allsub_HRFs_mean;

% p_diff_sr=nan(nROIs,nTRs_out);
for vv=1:nROIs
    for cc=1:nConds
        for tr=1:nTRs_out
            realvals = squeeze(vals(:,vv,cc,1:2,tr));
            % what is the sign-rank statistic for the real data?
            real_sr_stat(vv,cc,tr) = signrank_MMH(realvals(:,1),realvals(:,2));
            % determine before the parfor loop which conditions get randomly
            % swapped on each iteration (otherwise not deterministic)
            inds2swap = double(randn(nSubj,nPermIter)>0);
            inds2swap(inds2swap==0) = -1;
            parfor ii=1:nPermIter          
                % randomly permute the condition labels within subject
                randvals=realvals;
                randvals(inds2swap(:,ii)==-1,:) = randvals(inds2swap(:,ii)==-1,[2,1]);    
                % what is the sign-rank statistic for this randomly permuted data?
                rand_sr_stat(vv,cc,tr,ii) = signrank_MMH(randvals(:,1),randvals(:,2));
            end
        end
    end
end

% % compute a two-tailed p-value comparing the real stat to the random
% % distribution. Note that the <= and >= are inclusive, because any
% % iterations where real==null should count toward the null hypothesis. 
p_contraipsi_diff = 2*min(cat(4,mean(repmat(real_sr_stat,1,1,1,nPermIter)>=rand_sr_stat,4), ...
    mean(repmat(real_sr_stat,1,1,1,nPermIter)<=rand_sr_stat,4)),[],4);

array2table(squeeze(p_contraipsi_diff(vismotor_inds,1,1:15)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Pred_ContraIpsiDiff_TR',1:15))
array2table(squeeze(p_contraipsi_diff(vismotor_inds,1,16:30)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Pred_ContraIpsiDiff_TR',16:30))

array2table(squeeze(p_contraipsi_diff(vismotor_inds,2,1:15)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Rand_ContraIpsiDiff_TR',1:15))
array2table(squeeze(p_contraipsi_diff(vismotor_inds,2,16:30)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Rand_ContraIpsiDiff_TR',16:30))

%% compare contra-ipsi difference values between conditions
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 657567;
rng(rndseed,'twister')

real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);

for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,3,tr));

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

%% Now we make a figure for all subjects
sig_heights2 = [0.62, 0.655, 0.69];
if plotVisMotor
%     nplots_vis = ceil(sqrt(numel(vismotor_inds)));
    nplots_vis=5;
    for cc=1:nConds
        figure;hold all;
        for vi = 1:numel(vismotor_inds)
            subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
            ll=[];lx=0;lh=[];
            for rr=1:nResp
              
                vals = allsub_HRFs_mean(:,vismotor_inds(vi),cc,rr,:);
                if nSubj==1
                    meanvals =squeeze(vals);
                    semvals = squeeze(allsub_HRFs_sem(:,vismotor_inds(vi),cc,rr,:));
                else
                    meanvals = squeeze(nanmean(vals,1));
                    semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
                end

                lh=[lh,plot(tax,meanvals,'-','Color',col_resp(rr,:),'LineWidth',lw)];
                bandedError_MMH(tax, meanvals',semvals', col_resp(rr,:), 0.5);
%                     errorbar(tax, meanvals, semvals,'Color',col_resp(rr,:,bb),'LineWidth',lw);
                lx=lx+1;
                ll{lx}=sprintf('%s',respLabStrs{rr});
                
                
               
            end

%             for aa=[3]
            for aa=1:numel(alpha_vals)               
                inds2plot=p_contraipsi_diff(vismotor_inds(vi),cc,:)<alpha_vals(aa);
                plot(tax(inds2plot), repmat(sig_heights2(cc),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa));
                lh=[lh, plot(tax(1)-5, sig_heights2(cc),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
                ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
            end
                
            set(gca, 'FontSize', 12, 'XLim',[0 max(tax)],'YLim',ylims,'YTick',[0, 0.5])
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
        suptitle(sprintf('%s',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,600,1600]);
        saveas(gcf,fullfile(figpath,sprintf('ContraIpsi_%s.pdf',condLabStrs{cc})),'pdf');

    end
end

%% Now we make a figure for all subjects
sig_heights2 = [0.62, 0.655, 0.69];
subcolors = viridis(nSubj+1);
if plotVisMotorSS
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));

    for cc=1:nConds
        figure;hold all;
        for vi = 1:numel(vismotor_inds)
            subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
            ll=[];lx=0;lh=[];
            for ss=1:nSubj
              
                vals = squeeze(allsub_HRFs_mean(ss,vismotor_inds(vi),cc,3,:));
                
                lh=[lh,plot(tax,vals,'-','Color',subcolors(ss,:),'LineWidth',lw)];
                ll{ss}=sprintf('S%02d - contra-ipsi diff',sublist(ss));
                
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
        suptitle(sprintf('%s',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1800,1200]);
    end
end
%%
if plotVisMotor_Diff
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));
    figure;hold all;
      
    for vi = 1:numel(vismotor_inds)
        subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
        lh=[];lx=0;ll=[];
        for cc=1:nConds
          
            % take contra - ipsi difference
            vals = allsub_HRFs_mean(:,vismotor_inds(vi),cc,3,:);

            if nSubj==1
                meanvals =squeeze(vals);
                semvals = squeeze(allsub_HRFs_sem(:,vismotor_inds(vi),cc,3,:));
            else
                meanvals = squeeze(nanmean(vals,1));
                semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
            end

            lh=[lh,plot(tax,meanvals,'-','Color',col_conds(cc,:),'LineWidth',lw)];
            bandedError_MMH(tax, meanvals',semvals', col_conds(cc,:), 0.5);
%                 lh=[lh, errorbar(tax, meanvals, semvals,'Color',col_conds(cc,:,bb),'LineWidth',lw)];
            lx=lx+1;
            ll{lx} = sprintf('%s',condLabStrs{cc});

            for aa=1:numel(alpha_vals)
                inds2plot=p_contraipsi_diff(vismotor_inds(vi),cc,:)<alpha_vals(aa);
                plot(tax(inds2plot), repmat(sig_heights(cc),sum(inds2plot),1),'.','Color',col_conds(cc,:),'MarkerSize',alpha_ms(aa));
            end
            
        end
        
        for aa=1:numel(alpha_vals)
            inds2plot=p_cond_diff(vismotor_inds(vi),:)<alpha_vals(aa);
            plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa));
            lh=[lh, plot(tax(1)-5, sig_heights(nConds+1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
            ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
        end
                
            
        set(gca, 'FontSize', 12, 'XLim',[0 max(tax)],'YLim',difflims)
        h=plot(get(gca,'XLim'),[0,0],'--','Color',[0.7, 0.7, 0.7]);
        uistack(h,'bottom')
        for ee = 1:length(evts2plot)
            h=plot([evts2plot(ee),evts2plot(ee)],difflims,'-','Color',[0.8, 0.8, 0.8]);
            uistack(h,'bottom')
        end
        h=fill([evts2plot(1),evts2plot(2), evts2plot(2),evts2plot(1)], repelem(difflims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        h=fill([evts2plot(3),evts2plot(4), evts2plot(4),evts2plot(3)], repelem(difflims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
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
    suptitle('Contra finger - Ipsi finger')
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200]);
    
end
