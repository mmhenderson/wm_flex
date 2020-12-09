% script to plot the result of decoding analyses for oriSpin. 

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


plot_order = [12:14];  % motor areas
motor_names = ROI_names(plot_order);
plot_order_all = [plot_order];
motor_inds = find(ismember(plot_order_all,plot_order));
nROIs = length(plot_order_all);

plotMotorAcc=1;
plotMotorAccSS=1;
plotMotorDSS=1;

nVox2Use = 10000;
nPermIter=1000;
chance_val=0.5;
% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 1];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);
% cc=1;
alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,16,24];
alpha=alpha_vals(1);

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.93,0.96,0.99];
diff_col=[0.5, 0.5, 0.5];

fs=20;  % font size for all plots
ms=10;  % marker size for significance dots
%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
lw =1;


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);


acc_allsubs = nan(nSubj,nROIs,nConds,nTRs_out);
d_allsubs = nan(nSubj,nROIs,nConds,nTRs_out);

accrand_allsubs = nan(nSubj, nROIs, nConds, nTRs_out, nPermIter);


for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyActualResponse_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    acc_allsubs(ss,:,:,:) = squeeze(allacc(plot_order_all,:,:));
    d_allsubs(ss,:,:,:) = squeeze(alld(plot_order_all,:,:));
    
    accrand_allsubs(ss,:,:,:,:) = squeeze(allacc_rand(plot_order_all,:,:,:,:));
   
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
array2table(squeeze(p_sr(motor_inds,1,1:15)),...
    'RowNames',motor_names,'VariableNames',strseq('Pred_TR',1:15))
array2table(squeeze(p_sr(motor_inds,1,16:30)),...
    'RowNames',motor_names,'VariableNames',strseq('Pred_TR',16:30))
array2table(squeeze(p_sr(motor_inds,2,1:15)),...
    'RowNames',motor_names,'VariableNames',strseq('Rand_TR',1:15))
array2table(squeeze(p_sr(motor_inds,2,16:30)),...
    'RowNames',motor_names,'VariableNames',strseq('Rand_TR',16:30))

%% now doing pairwise condition comparisons - paired t-test.
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 877687;
rng(rndseed,'twister')

real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);
% p_diff_sr=nan(nROIs,nTRs_out);
for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,tr));
%         [p,h,stats]=signrank(realvals(:,1),realvals(:,2));
%         p_diff_sr(vv,tr) = p;

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
p_diff_sr = 2*min(cat(3,mean(repmat(real_sr_stat,1,1,nPermIter)>=rand_sr_stat,3), ...
    mean(repmat(real_sr_stat,1,1,nPermIter)<=rand_sr_stat,3)),[],3);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

array2table(squeeze(p_diff_sr(motor_inds,1:15)),...
    'RowNames',motor_names,'VariableNames',strseq('Diff_TR',1:15))
array2table(squeeze(p_diff_sr(motor_inds,16:30)),...
    'RowNames',motor_names,'VariableNames',strseq('Diff_TR',16:30))

%% Now we make a figure for all subjects

if plotMotorAcc
    nplots_vis = ceil(sqrt(numel(motor_inds)));
    figure();hold all;
    
    for vi = 1:numel(motor_inds)
       
        subplot(nplots_vis,ceil(numel(motor_inds)/nplots_vis),vi);hold all;
        
        lh=[];
        for cc = 1:nConds

            valsplot = acc_allsubs(:,motor_inds(vi),cc,:);
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
                inds2plot=p_sr(motor_inds(vi),cc,:)<alpha_vals(aa);
                plot(tax(inds2plot), repmat(sig_heights(cc),sum(inds2plot),1),'.','Color',col(cc,:),'MarkerSize',alpha_ms(aa));
            end
        end
        ll=condLabStrs;
        for aa=1:numel(alpha_vals)
            inds2plot=p_diff(motor_inds(vi),:)<alpha_vals(aa);
            plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))
            lh=[lh, plot(tax(1)-5, repmat(sig_heights(nConds+1),1,1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
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

        if vi==numel(motor_inds)
            legend(lh,ll,'FontSize', fs);
        end
       
        title(sprintf('%s', motor_names{vi}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
end
saveas(gcf,fullfile(figpath,'DecodeActualResp_allareas_overtime.pdf'),'pdf');

%%
if plotMotorAccSS
    sub_colors=viridis(nSubj+1);
    nplots_vis = ceil(sqrt(numel(motor_inds)));
    figure();hold all;
    lims=[0.3, 1];
    for vi = 1:numel(motor_inds)
       
        subplot(nplots_vis,ceil(numel(motor_inds)/nplots_vis),vi);hold all;

        cc=1;
        valsplot = acc_allsubs(:,motor_inds(vi),cc,:);
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

        if vi==numel(motor_inds)
            legend(lh,substrs,'FontSize', fs);
        end
       
        title(sprintf('%s', motor_names{vi}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
end
%%
if plotMotorDSS
    sub_colors=viridis(nSubj+1);
    nplots_vis = ceil(sqrt(numel(motor_inds)));
    figure();hold all;
    lims=[-1,3];
    for vi = 1:numel(motor_inds)
       
        subplot(nplots_vis,ceil(numel(motor_inds)/nplots_vis),vi);hold all;

        cc=1;
        valsplot = d_allsubs(:,motor_inds(vi),cc,:);
        lh=[];
        substrs=[];
        for ss=1:nSubj
            lh=[lh, plot(tax,squeeze(valsplot(ss,:)),'-','Color',sub_colors(ss,:),'LineWidth',lw)];
            substrs{ss} = sprintf('S%02d',sublist(ss));
        end

        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(lims);
        plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],lims,'-','Color',[0.8, 0.8, 0.8]);
        end
        
        if vi==1
            xlabel('Time(s)');
            ylabel('Dprime');
        end

        if vi==numel(motor_inds)
            legend(lh,substrs,'FontSize', fs);
        end
       
        title(sprintf('%s', motor_names{vi}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
end
