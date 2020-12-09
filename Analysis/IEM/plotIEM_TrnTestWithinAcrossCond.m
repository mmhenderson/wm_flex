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

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
       'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all','IPS0-3','IPS0-1','IPS2-3'};


plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
% plot_order1 = [1,22];
plot_order2 = [15:20];

visual_names = ROI_names(plot_order1);
% motor_names = ROI_names(plot_order2);
md_names = ROI_names(plot_order2);

plot_order_all = [plot_order1,plot_order2];
nROIs = length(plot_order_all);

vis_inds = find(ismember(plot_order_all,plot_order1));
% motor_inds = find(ismember(plot_order_all,plot_order2));
md_inds = find(ismember(plot_order_all,plot_order2));



nROIs = length(plot_order_all);


ylims_fid = [-0.1, 0.25];
col = plasma(5);
col = col(2:2:end-1,:);


condLabStrs = {'Trn/Test Predictable','Trn/Test Random','Trn Random/Test Predictable','Trn Predictable/Test Random'};
nConds = 2;

chance_val=0.5;

plotVisFids = 1;
plotVisFidSS=1;
plotVisFidAvg=1;
plotMDFid=1;
plotMDFidAvg=1;


diff_col=[0.5, 0.5, 0.5];
ms=10;  % marker size for significance dots
%% load results
nTrialsTotal = 400;

fid_allsubs = nan(nSubj,nROIs,nConds,2);    % last dim is trained on data from same or opp condition

condlabs_allsubs = nan(nSubj, nTrialsTotal);

for ss=1:length(sublist)
    
    substr = sprintf('S%02d',sublist(ss));    
    fn = fullfile(exp_path,'Analysis','IEM','IEM_results',sprintf('TrnTestWithinConds_%s.mat',substr));
    load(fn);
    assert(numel(allchanresp)==numel(ROI_names));
    allchanresp_within=allchanresp;
    
    fn = fullfile(exp_path,'Analysis','IEM','IEM_results',sprintf('TrnTestAcrossConds_%s.mat',substr));
    load(fn);
    assert(numel(allchanresp)==numel(ROI_names));
    allchanresp_across=allchanresp;
    
    condlabs = allchanresp_within(1).condLabs;
    assert(all(allchanresp_across(1).condLabs==condlabs));
    
    for vv=1:nROIs
        
       if plot_order_all(vv)>length(allchanresp) || isempty(allchanresp(plot_order_all(vv)).chan_resp_shift) 
           fprintf('skipping %s for S%s because no voxels\n', ROI_names{plot_order_all(vv)}, substr);
           continue
       end
       for cc = 1:nConds
           
           % take out trials from one condition at a time
           theserecs = allchanresp_within(plot_order_all(vv)).chan_resp_shift(condlabs==cc,:);           
           % get the fidelity
           angs = abs((1:360)-180);
           cos_vals = cosd(angs);
           fid_allsubs(ss,vv,cc,1) = mean(cos_vals.*mean(theserecs,1));
           
           % take out trials from one condition at a time
           theserecs = allchanresp_across(plot_order_all(vv)).chan_resp_shift(condlabs==cc,:);           
           % get the fidelity
           angs = abs((1:360)-180);
           cos_vals = cosd(angs);
           fid_allsubs(ss,vv,cc,2) = mean(cos_vals.*mean(theserecs,1));

       end
    end
    
end

assert(~any(isnan(fid_allsubs(:))))

%% More Stats
% for each decoding value, compute whether it's above chance across all
% subjects
% get a t-statistic: (mean-chance)/sem
vals = fid_allsubs;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
% tstat_real = (meanvals-chance_val)./semvals;
% 
% % same for the random values
% randvals = accrand_allsubs;
% meanvals_rand = squeeze(mean(randvals,1));
% semvals_rand = squeeze(std(randvals,[],1)./sqrt(nSubj));
% tstat_rand = (meanvals_rand-chance_val)./semvals_rand;
% 
% % finally get p-values based on how often real<random
% p = mean(permute(repmat(tstat_real,1,1,1,nPermIter),[1,2,4,3])<tstat_rand,3);
% is_sig = p<0.01;    % one tailed test
% 
% % ignore these results for areas we didn't actually test (have nans in the
% % random results marix)
% not_tested = squeeze(isnan(accrand_allsubs(1,:,1,1,2)));
% is_sig(not_tested,:,:,:) = 0;

% % now doing pairwise condition comparisons - paired t-test.
% numcores = 8;
% if isempty(gcp('nocreate'))
%     parpool(numcores);
% end
% rndseed = 345455;
% rng(rndseed,'twister');
% p_diff = nan(nROIs,1);
% for vv=1:nROIs
%     realvals = squeeze(vals(:,vv,:));
%     realdiff = diff(realvals,[],2);
%     realmeandiff = mean(realdiff);
%     inds2swap = double(randn(nSubj,nPermIter)>0);
%     inds2swap(inds2swap==0) = -1;
%     randmeandiffs = nan(nPermIter,1);
%     parfor ii=1:nPermIter    
%         % randomly swap which conditions are which, keeping subject labels
%         % same (this amounts to negating some of the differences)
%         randdiff = realdiff.*inds2swap(:,ii);
%         randmeandiffs(ii) = mean(randdiff);
%     end
%     % two-sided t-test here bc don't have obvious hypothesis
%     p_diff(vv) = 2*min([mean(realmeandiff<randmeandiffs), mean(realmeandiff>randmeandiffs)]);
% end
% diff_is_sig = p_diff<0.01;    

%% make a bar plot of acc - visual areas
col = gray(6);
col= col(1:4,:);
sub_colors = gray(nSubj+1);
if plotVisFids
    
    mean_plot = [meanvals(vis_inds,:,1),meanvals(vis_inds,:,2)];
    sem_plot = [semvals(vis_inds,:,1),semvals(vis_inds,:,2)];
%     sig_plot = [is_sig(vis_inds,:,1,1),is_sig(vis_inds,:,1,2)];
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,ylims_fid,visual_names,condLabStrs,...
        'Fidelity','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
    
    
end

%%
col = gray(3);
col= col(1:2,:);

if plotVisFidAvg
    
    vals = mean(fid_allsubs(:,vis_inds,:,:),3);
    mean_plot = squeeze(mean(vals,1));
    sem_plot = squeeze(std(vals,[],1)./sqrt(nSubj));
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,ylims_fid,visual_names,{'Within','Across'},...
        'Fidelity','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
    
    
end
%%
bw=0.50;
fs=14;
col = gray(6);
col= col(1:4,:);
sub_colors = plasma(nSubj+1);
if plotVisFidSS
   
    mean_plot = [meanvals(vis_inds,:,1),meanvals(vis_inds,:,2)];
    sem_plot = [semvals(vis_inds,:,1),semvals(vis_inds,:,2)];
%     sig_plot = [is_sig(vis_inds,:,1,1),is_sig(vis_inds,:,1,2)];
    
    

    set(groot,'DefaultLegendAutoUpdate','off');
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,ylims_fid,visual_names,condLabStrs,...
        'Fidelity','Spatial Memory Position',col)
    set(gca,'FontSize',fs);
    set(gcf,'Position',[800,800,1200,500])
    % get locations of bars w offsets
    c=get(gcf,'Children');b=get(c(end),'Children');
    bar_offsets = [b(end-3:end).XOffset];
%     set(b(3:10),'Visible','off')
    for vv=1:numel(vis_inds)
        for ss=1:nSubj
            subvals = [squeeze(fid_allsubs(ss,vis_inds(vv),:,1))', ...
                squeeze(fid_allsubs(ss,vis_inds(vv),:,2))'];
            h=plot(vv-bar_offsets,subvals,'.-','Color',sub_colors(ss,:),'LineWidth',1.5);
%             uistack(h,'bottom');
        end
    end
    b(end).BarWidth=bw;
    b(end-1).BarWidth=bw;
%     uistack(b(end),'top');
%     uistack(b(end-1),'top')
end


%% make a bar plot of acc - md areas
if plotMDFid
    
    mean_plot = [meanvals(md_inds,:,1),meanvals(md_inds,:,2)];
    sem_plot = [semvals(md_inds,:,1),semvals(md_inds,:,2)];
%     sig_plot = [is_sig(md_inds,:,1,1),is_sig(md_inds,:,1,2)];
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,ylims_fid,md_names,condLabStrs,...
        'Fidelity','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
end

%%
col = gray(3);
col= col(1:2,:);
sub_colors = gray(nSubj+1);
if plotMDFidAvg
    
    vals = mean(fid_allsubs(:,md_inds,:,:),3);
    mean_plot = squeeze(mean(vals,1));
    sem_plot = squeeze(std(vals,[],1)./sqrt(nSubj));
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,ylims_fid,md_names,{'Within','Across'},...
        'Fidelity','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
    
    
end