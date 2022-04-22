%% Plot accuracy of spatial decoder
% trained on spatial working memory mapping task and tested on each
% condition of main task. Time-resolved (TR by TR)
% Decoding analysis itself performed in TrnSWMLoc_TstWM_TRbyTR.m and saved 
% as mat file. This script loads that file, does all stats and plotting. 

% Makes Figure 2C, and Figure 2 - figure supplement 1.
%%
clear
close all;

sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;
filesepinds = find(curr_dir==filesep);
nDirsUp = 2;
exp_path = curr_dir(1:filesepinds(end-nDirsUp+1));
figpath = fullfile(exp_path,'figs');
addpath(fullfile(exp_path,'Analysis','stats_code'))
addpath(fullfile(exp_path,'Analysis','plotting_utils'));

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

% Indices into "ROI_names" corresponding to visual ROIs and motor ROIs
plot_order = [1:5,10,11,6:9,12:14];  
vis_names = ROI_names(plot_order);
plot_order_all = [plot_order];
vis_inds = find(ismember(plot_order_all,plot_order));
nROIs = length(plot_order_all);

plotVisAcc=1;   % plot subject-averaged decoding acc w error bars?
plotVisAccSS=1;  % plot single-subject decoding performance for each cond over time?

nVox2Use = 10000; 
nPermIter=1000;
class_str = 'normEucDist';

% plotting and stats info
acclims = [0.4, 1];
col = [3, 70, 124; 110, 172, 229]./255;

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);
chance_val=0.5;

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.90,0.94,0.98];
diff_col=[0.5, 0.5, 0.5];

fs=12;  % font size for all plots
ms=10;  % marker size for significance dots
%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
lw =1;


condLabStrs = {'Informative','Uninformative'};
nConds = length(condLabStrs);


acc_allsubs = nan(nSubj,nROIs,nConds,nTRs_out);
d_allsubs = nan(nSubj,nROIs,nConds,nTRs_out);

accrand_allsubs = nan(nSubj, nROIs, nConds, nTRs_out, nPermIter);


for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_TRbyTR_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    % averaging decoding performance over the four decoding schemes (0 vs
    % 180, 45 versus 225, etc)
    acc_allsubs(ss,:,:,:) = squeeze(mean(allacc(plot_order_all,:,:,:),3));
    d_allsubs(ss,:,:,:) = squeeze(mean(alld(plot_order_all,:,:,:),3));
    
    accrand_allsubs(ss,:,:,:,:) = squeeze(mean(allacc_rand(plot_order_all,:,:,:,:,:),3));
   
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

% compute some basic stats here
vals = acc_allsubs;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
randvals = accrand_allsubs;

% which areas did we do significance test on? not all, because it was slow
inds2test = ~isnan(squeeze(accrand_allsubs(1,:,1,1,1)));

%% Wilcoxon signed rank test
% for each permutation iteration, use this test to compare real data for all subj to
% shuffled data for all subj.
stat_iters_sr = nan(nROIs, nConds, nTRs_out, nPermIter); 

for vv=1:nROIs
    for cc=1:nConds
        for tr = 1:nTRs_out
            x = vals(:,vv,cc,tr);
        
            for ii=1:nPermIter
                y = randvals(:,vv,cc,tr,ii);
                % compare the real values against the null, for this iteration.
                % w>0 means real>null, w<0 means real<null, w=0 means equal
                stat_iters_sr(vv,cc,tr,ii) = signrank_MMH(x,y);

            end
        end
    end
end

% final p value is the proportion of iterations where null was at least as
% large as the real (e.g. the test stat was 0 or negative)
p_sr = mean(stat_iters_sr<=0, 4);
% for any areas we didn't actually test (nans in the shuffled accuracy), 
% un-mark as significant
p_sr(~inds2test,:,:) = 100;

% print out p values for each condition
array2table(squeeze(p_sr(vis_inds,1,1:15)),...
    'RowNames',vis_names,'VariableNames',strseq('Pred_TR',1:15))
array2table(squeeze(p_sr(vis_inds,1,16:30)),...
    'RowNames',vis_names,'VariableNames',strseq('Pred_TR',16:30))
array2table(squeeze(p_sr(vis_inds,2,1:15)),...
    'RowNames',vis_names,'VariableNames',strseq('Rand_TR',1:15))
array2table(squeeze(p_sr(vis_inds,2,16:30)),...
    'RowNames',vis_names,'VariableNames',strseq('Rand_TR',16:30))

%% pairwise condition comparisons within each timept
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 745675;
rng(rndseed,'twister')

real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);

for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,tr));

        % what is the sign-rank statistic for the real data?
        real_sr_stat(vv,tr) = signrank_MMH(realvals(:,1),realvals(:,2));

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

% compute a two-tailed p-value comparing the real stat to the random
% distribution. Note that the <= and >= are inclusive, because any
% iterations where real==null should count toward the null hypothesis. 
p_diff_sr = 2*min(cat(3,mean(repmat(real_sr_stat,1,1,nPermIter)>=rand_sr_stat,3), ...
    mean(repmat(real_sr_stat,1,1,nPermIter)<=rand_sr_stat,3)),[],3);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

array2table(squeeze(p_diff_sr(vis_inds,1:15)),...
    'RowNames',vis_names,'VariableNames',strseq('Diff_TR',1:15))
array2table(squeeze(p_diff_sr(vis_inds,16:30)),...
    'RowNames',vis_names,'VariableNames',strseq('Diff_TR',16:30))

%% Now we make a figure for all subjects

if plotVisAcc

    nplots_vis = ceil(sqrt(numel(vis_inds)));
    figure();hold all;
    
    for vi = 1:numel(vis_inds)
       
        subplot(nplots_vis,ceil(numel(vis_inds)/nplots_vis),vi);hold all;
        
        lh=[];
        for cc = 1:nConds

            valsplot = acc_allsubs(:,vis_inds(vi),cc,:);
            if nSubj==1
                meanvals =squeeze(valsplot)';
                semvals = [];
            else
                meanvals = squeeze(mean(valsplot,1))';
                semvals = squeeze(std(valsplot,[],1))'./sqrt(nSubj);
            end
            lh=[lh,plot(tax,meanvals,'-','Color',col(cc,:),'LineWidth',lw)];
            bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);

            for aa=1:numel(alpha_vals)
                inds2plot=p_sr(vis_inds(vi),cc,:)<alpha_vals(aa);
                plot(tax(inds2plot), repmat(sig_heights(cc),sum(inds2plot),1),'.','Color',col(cc,:),'MarkerSize',alpha_ms(aa));
            end
        end
        ll=condLabStrs;
        for aa=1:numel(alpha_vals)
            inds2plot=p_diff(vis_inds(vi),:)<alpha_vals(aa);
            plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))
            lh=[lh, plot(tax(1)-5, repmat(sig_heights(nConds+1),1,1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
            ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
        end
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(acclims);
        h=plot(get(gca,'XLim'),[chance_val,chance_val],'--','Color',[0.7, 0.7, 0.7]);
        uistack(h,'bottom')
        for ee = 1:length(evts2plot)
            h=plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
            uistack(h,'bottom')
        end
        h=fill([evts2plot(1),evts2plot(2), evts2plot(2),evts2plot(1)], repelem(acclims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        h=fill([evts2plot(3),evts2plot(4), evts2plot(4),evts2plot(3)], repelem(acclims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        
        if vi==1
            xlabel('Time(s)');
            ylabel('Accuracy');
        end

        if vi==numel(vis_inds)
            legend(lh,ll,'FontSize', fs);
        end
       
        title(sprintf('%s', vis_names{vi}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200]);
end
saveas(gcf,fullfile(figpath,'TrainSWM_TestWM_allareas_overtime.pdf'),'pdf');

%%
if plotVisAccSS
    for cc=1:nConds
        sub_colors=viridis(nSubj+1);
        nplots_vis = ceil(sqrt(numel(vis_inds)));
        figure();hold all;
        lims=[0.3, 1];
        for vi = 1:numel(vis_inds)

            subplot(nplots_vis,ceil(numel(vis_inds)/nplots_vis),vi);hold all;


            valsplot = acc_allsubs(:,vis_inds(vi),cc,:);
            lh=[];
            substrs=[];
            for ss=1:nSubj
                lh=[lh, plot(tax,squeeze(valsplot(ss,:)),'-','Color',sub_colors(ss,:),'LineWidth',lw)];
                substrs{ss} = sprintf('S%02d',sublist(ss));
            end

            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(lims);
            plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],lims,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Accuracy');
            end

            if vi==numel(vis_inds)
                legend(lh,substrs,'FontSize', fs);
            end

            title(sprintf('%s', vis_names{vi}));

        end
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1800,1200]);
        suptitle(condLabStrs{cc})
    end
end
