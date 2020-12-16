%% Plot average signal in each ROI over time. 
% similar to deconvolution but instead of deconvolving, just stack epoched
% trials and take average. 

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

plotVisMotor = 1;
plotVisMotorSS = 1;

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
nROIs = length(plot_order_all);

lw=1;

% tr2baseline=0;

condLabStrs = {'Predictable','Random'};

% length of each run type
nTRs = 583-16;
trDur = 0.80;

% deconvolution param
nTRs_out = 30; % # of TRs post-stimulus to model in the decon analysis
% define timing axis
tax = trDur*(0:nTRs_out-1); 

nTrialsPerRun = 20;
nConds = 2;

allsub_HRFs_mean = nan(nSubj, length(plot_order_all), nConds, nTRs_out);
allsub_HRFs_sem = nan(nSubj, length(plot_order_all), nConds, nTRs_out);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

ylims=[-0.5, 0.5];
sig_heights = [0.42,0.455,0.49];

col_conds = [125, 93, 175; 15, 127, 98]./255;

diff_col=[0.5, 0.5, 0.5];
%%
for ss=1:length(sublist)

    
    substr = sprintf('S%02d',sublist(ss));
   
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    
    %
    for vi = 1:nROIs % for all visual areas I want to look at
        vv = plot_order_all(vi);
        
        % trials by TRs by voxels
        dat_by_TR = mainSig(vv).dat_by_TR;
        % baseline to zero at first TR
%         dat_by_TR = dat_by_TR - repmat(dat_by_TR(:,1,:),1,size(dat_by_TR,2),1);
        
        condLabs = mainSig(vv).condLabs;
        nVox = size(dat_by_TR,3);
        
        % take an average over voxels and save it for later
        for cc = 1:nConds
           mean_over_trials= squeeze(mean(dat_by_TR(condLabs==cc,:,:),1));
           allsub_HRFs_mean(ss,vi,cc,:) = mean(mean_over_trials,2);
           allsub_HRFs_sem(ss,vi,cc,:) = std(mean_over_trials,[],2)./sqrt(size(mean_over_trials,2));
        end
       
        if nVox==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end

        fprintf('processing %s area %s, %d voxels\n', substr, ROI_names{vv}, nVox);
        
      
 
    end

end


%% compare between conditions
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 456556;
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
    suptitle('Mean Z-scored BOLD, all trials each condition')
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200]);
    
    
end
%%
subcolors = viridis(nSubj+1);
sslims =  [-0.7, 0.7];
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
        suptitle(sprintf('%s - mean Z-scored BOLD over all trials',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1800,1200]);
    end
end
