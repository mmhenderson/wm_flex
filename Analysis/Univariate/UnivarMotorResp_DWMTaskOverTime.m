%% Plot signal related to contralateral/ipsilateral finger presses
% during the digit working memory (DWM) task, delayed button pressing task.
% Get the average univariate signal in each hemisphere ROI on trials
% where either the left or right finger was to be pressed at end of delay. 
% See when/where univariate motor signals emerge.
% Note this task is not included in our paper.
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

% names of the ROIs - note the motor areas have been changed to nicer names
% here for plotting.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

plotVisMotor=1;
plotVisMotorSS=1;
plotVisMotorDiff=1;

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));

nROIs = length(plot_order_all);
lw=1;

condLabStrs = {'DWM Loc Task'};
respLabStrs = {'Contra finger','Ipsi finger'};

% length of each run type
nTRs = 452-16;
trDur = 0.80;

% deconvolution param
nTRs_out = 25; % # of TRs post-stimulus to model in the decon analysis
% define timing axis
tax = trDur*(0:nTRs_out-1); 

nTrialsPerRun = 20;
nConds = 1;
nResp=2;

% the nResp dimension goes [contra, ipsi, contra - ipsi]
allsub_HRFs_mean = zeros(nSubj, length(plot_order_all), nConds, nResp+1, nTRs_out);
allsub_HRFs_sem = zeros(nSubj, length(plot_order_all), nConds, nResp+1, nTRs_out);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);

% events to plot as vertical lines (boundary onset/offsets)
evts2plot = [0,1,13,14];

ylims = [-0.5, 0.7];
difflims = [-0.5 0.5];

sig_heights = [0.42,0.455,0.49];
diff_col=[0.5, 0.5, 0.5];

col_resp = plasma(3);
col_resp = col_resp(1:2,:);

col_conds = viridis(4);
col_conds=col_conds(1,:);

%%
for ss=1:length(sublist)

    
    substr = sprintf('S%02d',sublist(ss));
   
    fn2load = fullfile(exp_path,'Samples',sprintf('DWMLocSignalByTrial_sepHemis_%s.mat',substr));
    load(fn2load);
    
    %%
    for vi = 1:nROIs % for all visual areas I want to look at
        
        vv = plot_order_all(vi);



        % Pull out data from both sides, separately
        % trials by TRs by voxels
        dat_by_TR_lh = locSig(vv,1).dat_by_TR;
        dat_by_TR_rh = locSig(vv,2).dat_by_TR;
        if numel(dat_by_TR_lh)==0 && numel(dat_by_TR_rh)==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end
       
        respLabs = locSig(vv,1).ExpDigit;
        condLabs = ones(size(respLabs));
        % take an average over voxels and save it for later
        for cc = 1:nConds
            for rr=1:nResp
                % rr=1 for contra 
                    % want data from lh in brain, when right finger (respLabs==2) was used)
                    % want data from rh in brain, when left finger (respLabs==1) was used)
                % rr=2 for ipsi 
                    % want data from lh in brain, when left finger (respLabs==1) was used)
                    % want data from rh in brain, when right finger (respLabs==2) was used)
                if numel(dat_by_TR_lh)>0
                    dat1 = squeeze(dat_by_TR_lh(condLabs==cc & respLabs==(3-rr),:,:));
                else
                    dat1=[];
                end
                if numel(dat_by_TR_rh)>0
                    dat2 = squeeze(dat_by_TR_rh(condLabs==cc & respLabs==rr,:,:));
                else
                    dat2=[];
                end
                dat = cat(3,dat1,dat2);
                
                mean_over_trials = squeeze(mean(dat,1));
                if rr==1
                    contradat_mean_over_trials = mean_over_trials;
                elseif rr==2
                    ipsidat_mean_over_trials = mean_over_trials;
                end
                allsub_HRFs_mean(ss,vi,cc,rr,:) = mean(mean_over_trials,2);
                allsub_HRFs_sem(ss,vi,cc,rr,:) = std(mean_over_trials,[],2)./sqrt(size(mean_over_trials,2));
            end
            
            % now take the contra-ipsi difference
            diffdat = contradat_mean_over_trials-ipsidat_mean_over_trials;
            
            allsub_HRFs_mean(ss,vi,cc,3,:) = mean(diffdat,2);
            allsub_HRFs_sem(ss,vi,cc,3,:) = std(diffdat,[],2)./sqrt(size(diffdat,2));
            
        end
        fprintf('processing %s area %s\n', substr, ROI_names{vv});

    end
    

end


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
    'RowNames',vismotor_names,'VariableNames',strseq('DWM_ContraIpsiDiff_TR',1:15))
array2table(squeeze(p_contraipsi_diff(vismotor_inds,1,16:25)),...
    'RowNames',vismotor_names,'VariableNames',strseq('DWM_ContraIpsiDiff_TR',16:25))

%% Now we make a figure for all subjects
sig_heights2 = [0.62, 0.655, 0.69];
if plotVisMotor
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));

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

            for aa=1:numel(alpha_vals)
                inds2plot=p_contraipsi_diff(vismotor_inds(vi),cc,:)<alpha_vals(aa);
                plot(tax(inds2plot), repmat(sig_heights2(cc),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa));
                lh=[lh, plot(tax(1)-5, sig_heights2(cc),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
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
if plotVisMotorDiff
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
%             inds2plot=p_cond_diff(vismotor_inds(vi),:)<alpha_vals(aa);
%             plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa));
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
