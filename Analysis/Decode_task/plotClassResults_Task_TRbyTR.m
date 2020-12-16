%% Plot accuracy of task decoder
% trained on main task data, within each TR (time resolved)
% Decoding analysis itself performed in Classify_Task_TRbyTR.m and saved as mat
% file. This script loads that file, does all stats and plotting. 
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
addpath(fullfile(exp_path,'Analysis','stats_code','bayesian-prevalence','matlab'))

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

nVox2Use = 10000;
nPermIter=1000;
chance_val=0.5;

class_str = 'normEucDist';

% plotting and stats info
acclims = [0.4, 1];
dprimelims = [-0.2, 1.4];
col = viridis(4);
col = col(2,:);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,16,24];
alpha=alpha_vals(1);

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.93,0.96,0.99];

fs=20;  % font size for all plots
ms=10;  % marker size for significance dots
%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
lw =1;

acc_allsubs = nan(nSubj,nROIs,nTRs_out);
d_allsubs = nan(nSubj,nROIs,nTRs_out);

accrand_allsubs = nan(nSubj, nROIs, nTRs_out, nPermIter);


for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyTask_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    acc_allsubs(ss,:,:) = squeeze(allacc(plot_order_all,:));
    d_allsubs(ss,:,:) = squeeze(alld(plot_order_all,:));
    
    accrand_allsubs(ss,:,:,:) = squeeze(allacc_rand(plot_order_all,:,:));
   
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

% compute some basic stats here
vals = acc_allsubs;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
randvals = accrand_allsubs;

% which areas did we do significance test on? not all, because it was slow
inds2test = ~isnan(squeeze(accrand_allsubs(1,:,1,1)));

%% Wilcoxon signed rank test
% for each permutation iteration, use this test to compare real data for all subj to
% shuffled data for all subj.
stat_iters_sr = nan(nROIs, nTRs_out, nPermIter); 

for vv=1:nROIs
 
    for tr = 1:nTRs_out
        x = vals(:,vv,tr);

        for ii=1:nPermIter
            y = randvals(:,vv,tr,ii);
            % compare the real values against the null, for this iteration.
            % w>0 means real>null, w<0 means real<null, w=0 means equal
            stat_iters_sr(vv,tr,ii) = signrank_MMH(x,y);

        end
    end

end

% final p value is the proportion of iterations where null was at least as
% large as the real (e.g. the test stat was 0 or negative)
p_sr = mean(stat_iters_sr<=0, 3);
% for any areas we didn't actually test (nans in the shuffled accuracy), 
% un-mark as significant
p_sr(~inds2test,:) = 100;

% print out p values for each condition
array2table(squeeze(p_sr(vis_inds,1:15)),...
    'RowNames',vis_names,'VariableNames',strseq('Task_TR',1:15))
array2table(squeeze(p_sr(vis_inds,16:30)),...
    'RowNames',vis_names,'VariableNames',strseq('Task_TR',16:30))

%% Now we make a figure for all subjects
cc=1;
if plotVisAcc
    nplots_vis = ceil(sqrt(numel(vis_inds)));
    figure();hold all;
    
    for vi = 1:numel(vis_inds)
       
        subplot(nplots_vis,ceil(numel(vis_inds)/nplots_vis),vi);hold all;
        
        lh=[];
       
        valsplot = acc_allsubs(:,vis_inds(vi),:);
        if nSubj==1
            meanvals =squeeze(valsplot)';
            semvals = [];
        else
            meanvals = squeeze(mean(valsplot,1))';
            semvals = squeeze(std(valsplot,[],1))'./sqrt(nSubj);
        end
        lh=[lh,plot(tax,meanvals,'-','Color',col(cc,:),'LineWidth',lw)];
        bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);

        ll={'Task'};
        for aa=1:numel(alpha_vals)
            inds2plot=p_sr(vis_inds(vi),:)<alpha_vals(aa);
            h=plot(tax(inds2plot), repmat(sig_heights(cc),sum(inds2plot),1),'.','Color',col(cc,:),'MarkerSize',alpha_ms(aa));
            lh=[lh,h];
            ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
        end
     
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(acclims);
        plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
        end
        
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
