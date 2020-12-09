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


plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs and motor ROIs
plot_order2 = [15:20];  % MD areas (not super interested in these)

vismotor_names = ROI_names(plot_order1);
md_names = ROI_names(plot_order2);

plot_order_all = [plot_order1,plot_order2];
nROIs = length(plot_order_all);

vismotor_inds = find(ismember(plot_order_all,plot_order1));
md_inds = find(ismember(plot_order_all,plot_order2));


nVox2Use = 10000;
nPermIter=1000;
chance_val=0.5;
% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

alpha_vals=[0.05, 0.01, 0.001];
alpha_ms = [8,16,24];
alpha = alpha_vals(1);

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);

diff_col=[0.5, 0.5, 0.5];

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

plotVisMotorAcc = 1;
plotPrevalence=1;
plotMDAcc=1;
%% load results

acc_allsubs = nan(nSubj,nROIs,nConds);
accrand_allsubs = nan(nSubj,nROIs,nConds,nPermIter);
d_allsubs = nan(nSubj,nROIs,nConds);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_CorrectOnly_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
%     fn2load = fullfile(save_dir,sprintf('ClassifyResponse_leavePairOut_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    
    assert(size(allacc,1)==numel(ROI_names))
    acc_allsubs(ss,:,:) = squeeze(allacc(plot_order_all,:));
    d_allsubs(ss,:,:) = squeeze(alld(plot_order_all,:));
    
    accrand_allsubs(ss,:,:,:) = squeeze(allacc_rand(plot_order_all,:,:));
   
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

% get some basic stats to use for the plots and tests below
vals = acc_allsubs;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
randvals = accrand_allsubs;


%% 2-way RM anova on decoding values
% using shuffling to compute significance of each effect
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 756765;
[p_vals, ranova_table, iter] = get_f_dist(acc_allsubs(:,vismotor_inds,:), nPermIter, rndseed, 0);

% print results of the shuffling test
f_vals = ranova_table{[3,5,7],[4]}';
df = ranova_table{[3,5,7],[2]}';
array2table([f_vals; df; p_vals],'RowNames',{'f','df','pval'},'VariableNames',{'ROI','Condition','interaction'})

%% Wilcoxon signed rank test
% for each permutation iteration, use this test to compare real data for all subj to
% shuffled data for all subj.
stat_iters_sr = nan(nROIs, nConds, nPermIter); 

for vv=1:nROIs
    for cc=1:nConds
        x = vals(:,vv,cc);
        for ii=1:nPermIter
            y = randvals(:,vv,cc,ii);
            % compare the real values against the null, for this iteration.
            % w>0 means real>null, w<0 means real<null, w=0 means equal
            stat_iters_sr(vv,cc,ii) = signrank_MMH(x,y);

        end
    end
end

% final p value is the proportion of iterations where null was at least as
% large as the real (e.g. the test stat was 0 or negative)
p_sr = mean(stat_iters_sr<=0, 3);

is_sig=p_sr<alpha;

% print out how many which areas and conditions are significant across all
% subs
array2table([squeeze(is_sig(vismotor_inds,1)), p_sr(vismotor_inds,1), ...
    squeeze(is_sig(vismotor_inds,2)), p_sr(vismotor_inds,2)],...
    'RowNames',vismotor_names,'VariableNames',{'Pred','pval_pred','Rand','pval_rand'})

%% now doing pairwise condition comparisons - paired t-test.
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 123344;
rng(rndseed,'twister')
% 
real_sr_stat = nan(nROIs,1);
rand_sr_stat = nan(nROIs, nPermIter);
% p_diff_sr=nan(nROIs,1);
for vv=1:nROIs
    realvals = squeeze(vals(:,vv,:));
%     [p,h,stats]=signrank(realvals(:,1),realvals(:,2));
%     p_diff_sr(vv) = p;
    % what is the sign-rank statistic for the real data?
    real_sr_stat(vv) = signrank_MMH(realvals(:,1),realvals(:,2));
%     
    % determine before the parfor loop which conditions get randomly
    % swapped on each iteration (otherwise not deterministic)
    inds2swap = double(randn(nSubj,nPermIter)>0);
    inds2swap(inds2swap==0) = -1;

    parfor ii=1:nPermIter          
        
        % randomly permute the condition labels within subject
        randvals=realvals;
        randvals(inds2swap(:,ii)==-1,:) = randvals(inds2swap(:,ii)==-1,[2,1]);    
        % what is the sign-rank statistic for this randomly permuted data?
        rand_sr_stat(vv,ii) = signrank_MMH(randvals(:,1),randvals(:,2));

    end
end

% % compute a two-tailed p-value comparing the real stat to the random
% % distribution. Note that the <= and >= are inclusive, because any
% % iterations where real==null should count toward the null hypothesis. 
p_diff_sr = 2*min([mean(repmat(real_sr_stat,1,nPermIter)>=rand_sr_stat,2), ...
    mean(repmat(real_sr_stat,1,nPermIter)<=rand_sr_stat,2)],[],2);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

% print out which areas show a significant condition effect across all subj
array2table([diff_is_sig(vismotor_inds), p_diff(vismotor_inds)],'RowNames',vismotor_names,'VariableNames',{'cond_diff','p'})

%% compute individual subject significance of decoding
vals = acc_allsubs;
randvals = accrand_allsubs;
% finally get p-values based on how often real<random
p_ss = mean(repmat(vals,1,1,1,nPermIter)<randvals,4);
is_sig_ss = p_ss<alpha;    % one tailed test

% print out how many subjects were individually significant for each area
array2table(squeeze(sum(is_sig_ss(:,vismotor_inds,:),1)),'RowNames',vismotor_names,'VariableNames',condLabStrs)

% print out how many subjects individually showed condition difference in
% the direction of the main effect (pred>random
array2table(squeeze(sum(vals(:,vismotor_inds,1)>vals(:,vismotor_inds,2),1))','RowNames',vismotor_names,'VariableNames',{'pred_gr_rand'})


%% bayesian prevalance of significant effects in the population
n=nSubj;
a=alpha;
b=1;
prev=zeros(nROIs,nConds);
prev_lb=zeros(nROIs,nConds);
prev_hpdi95= zeros(nROIs, nConds, 2);
prev_hpdi50= zeros(nROIs, nConds, 2);
for vv=1:nROIs
    for cc=1:nConds
        k=sum(is_sig_ss(:,vv,cc));        
        prev(vv,cc) = bayesprev_map(k,n,a,b);
        post = bayesprev_posterior(0:0.05:1, k, n, a, b);
        prev_hpdi95(vv,cc,:) = bayesprev_hpdi(0.95, k, n, a, b);
        prev_hpdi50(vv,cc,:) = bayesprev_hpdi(0.50, k, n, a, b);
        prev_lb(vv,cc) = bayesprev_bound(0.95,k,n,a,b);
        % note this func only does lower bound, not upper
    end
end

% compare difference in prevalence between conditions
prev_bw = zeros(nROIs,1);
prev_hpdi95_bw= zeros(nROIs, 2);
prev_hpdi50_bw= zeros(nROIs, 2);
for vv=1:nROIs
    k11=sum(is_sig_ss(:,vv,1) & is_sig_ss(:,vv,2));   
    k10=sum(is_sig_ss(:,vv,1) & ~is_sig_ss(:,vv,2));   
    k01=sum(~is_sig_ss(:,vv,1) & is_sig_ss(:,vv,2));   
    [map, post_x, post_p, hpi, probGT, logoddsGT, samples] = bayesprev_diff_within(k11, k10, k01, n, 0.95, a, b, 10000);
    prev_bw(vv) = map;
    prev_hpdi95_bw(vv,:) = hpi;
    [map, post_x, post_p, hpi, probGT, logoddsGT, samples] = bayesprev_diff_within(k11, k10, k01, n, 0.50, a, b, 10000);
    prev_hpdi50_bw(vv,:) = hpi;
end

%%
if plotPrevalence
    figure;hold all;
    lw96= 4;
    lw50 = 8;
    medvals=prev;
    xpos = 1:numel(vismotor_inds);
    xoffset=[-0.15, 0.15];
    ylim([-0.2, 1.2])
    xlim([0, numel(vismotor_inds)+1]);
    plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8])
    lh=[];
    for vv=1:numel(vismotor_inds)
        for cc=1:nConds
            plot(repmat([xpos(vv)+xoffset(cc)],2,1), [prev_hpdi95(vismotor_inds(vv),cc,1),prev_hpdi95(vismotor_inds(vv),cc,2)],'-','LineWidth',lw96,'Color',col(cc,:));
            h=plot(repmat([xpos(vv)+xoffset(cc)],2,1), [prev_hpdi50(vismotor_inds(vv),cc,1),prev_hpdi50(vismotor_inds(vv),cc,2)],'-','LineWidth',lw50,'Color',col(cc,:));
            plot(xpos(vv)+xoffset(cc), medvals(vismotor_inds(vv),cc),'.','MarkerSize',12,'Color','k');
            if vv==1
                lh=[lh, h];
            end
        end
    end
    
    set(gcf,'Color','w');
    set(gca,'XTick',1:numel(vismotor_inds),'XTickLabels',vismotor_names,'XTickLabelRotation',90)
    ylabel('Population prevalence')
    title('Prevalence of significant response decoding in each condition')
    legend(lh,condLabStrs)
    
    % plot difference in prevalence
    figure;hold all;
    lw96= 4;
    lw50 = 8;
    medvals=prev_bw;
    xpos = 1:numel(vismotor_inds);
    ylim([-0.5, 1])
    xlim([0, numel(vismotor_inds)+1]);
    plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8])
    for vv=1:numel(vismotor_inds)
        plot(repmat([xpos(vv)],2,1), [prev_hpdi95_bw(vismotor_inds(vv),1),prev_hpdi95_bw(vismotor_inds(vv),2)],'-','LineWidth',lw96,'Color',diff_col);
        h=plot(repmat([xpos(vv)],2,1), [prev_hpdi50_bw(vismotor_inds(vv),1),prev_hpdi50_bw(vismotor_inds(vv),2)],'-','LineWidth',lw50,'Color',diff_col);
        plot(xpos(vv)+xoffset(cc), medvals(vismotor_inds(vv)),'.','MarkerSize',12,'Color','k');

    end
    
    set(gcf,'Color','w');
    set(gca,'XTick',1:numel(vismotor_inds),'XTickLabels',vismotor_names,'XTickLabelRotation',90)
    ylabel('Population prevalence')
    title('Estimated prevalence difference (rand-pred)')
    legend(lh,condLabStrs)
end

 %% plot with single subjects
bw=0.50;
fs=14;
if plotVisMotorAcc
   
    meanVals=meanvals(vismotor_inds,:);
    seVals=semvals(vismotor_inds,:);
    
    sub_colors = gray(nSubj+1);
    set(groot,'DefaultLegendAutoUpdate','off');
    fh = figure();hold on;
    % first make the actual bar plot
    b = bar(gca,meanVals);
    lh=[b(1),b(2)];
    
    % have to set this to "modal", otherwise it fails to get the XOffset
    % property.
    set(fh, 'WindowStyle','modal','WindowState','minimized')
    bar_offset = [b.XOffset];
    barPos = repmat((1:size(meanVals,1))', 1, length(bar_offset)) + repmat(bar_offset, size(meanVals,1), 1);
    for cc=1:nConds
        b(cc).FaceColor = col(cc,:);
        b(cc).EdgeColor = col(cc,:);
        errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    end

    set(gca,'XTick', 1:numel(vismotor_inds))
    set(gca,'XTickLabel', vismotor_names,'XTickLabelRotation',90);
    ylabel('Accuracy')
    set(gca,'YLim',acclims)
    set(gca,'XLim',[0,numel(vismotor_inds)+1])
    if chance_val~=0
        line([0,numel(vismotor_inds)+1],[chance_val,chance_val],'Color','k');
    end
    set(gca,'FontSize',fs);
    set(gcf,'Position',[800,800,1200,500]);
    % get locations of bars w offsets
    c=get(gcf,'Children');b=get(c(end),'Children');
   
    verspacerbig = range(acclims)/50;
    horspacer = abs(diff(bar_offset))/2;
%     
    for vv=1:numel(vismotor_inds)
        % add individual subjects
        for ss=1:nSubj
            subvals = squeeze(acc_allsubs(ss,vismotor_inds(vv),:));
            h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(5,:),'LineWidth',1.5);
            uistack(h,'bottom');
        end
        % add significance of individual areas/conditions
        for cc=1:nConds
            for aa=1:numel(alpha_vals)
                if p_sr(vv,cc)<alpha_vals(aa)
                    % smaller dots get over-drawn with larger dots
                    plot(vv+bar_offset(cc), meanVals(vv,cc)+seVals(vv,cc)+verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa))
                end
            end
        end
        % add significance of condition differences
        for aa=1:numel(alpha_vals)
            if p_diff(vv)<alpha_vals(aa)
                [mx,maxind] = max(meanVals(vv,:));
                % smaller dots get over-drawn with larger dots
                plot(vv+bar_offset, repmat(meanVals(vv,maxind)+seVals(vv,maxind)+2*verspacerbig,2,1),'-','Color','k','LineWidth',1)
                plot(vv, meanVals(vv,maxind)+seVals(vv,maxind)+3*verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa));
                
            end
            if vv==1
                lh=[lh,plot(-1, meanVals(vv,1)+seVals(vv,1)+3*verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa))];
            end
        end
    end
    b(end).BarWidth=bw;
    b(end-1).BarWidth=bw;
    leg=legend(lh,{'Predictable','Random','p<0.05','0<0.01','p<0.001'},'Location','EastOutside');
    set(gcf,'color','white')
    set(gcf, 'WindowStyle','normal','WindowState','normal')

end

saveas(gcf,fullfile(figpath,'DecodeCorrectResp_allareas.pdf'),'pdf');

%% make a bar plot of acc - md areas
if plotMDAcc
    
    plot_barsAndStars(meanvals(md_inds,:),semvals(md_inds,:),is_sig(md_inds,:),...
        diff_is_sig(md_inds),chance_val,acclims,md_names,condLabStrs,...
        'Accuracy','Response (expected)',col);
end
