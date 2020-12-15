%% Plot accuracy of spatial decoder
% trained and tested within OR across conditions of main task. 
% want to see if decoding is impaired when across versus within conditions.
% Decoding analysis itself performed in TrnWithinCond_leavePairOut.m and 
% TrnAcrossCond_leavePairOut.m, and saved as mat file. 
% This script loads those files, makes plots. 
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

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
       'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all','IPS0-3','IPS0-1','IPS2-3'};


% Indices into "ROI_names" corresponding to visual ROIs and motor ROIs
% reordering them a little for logical x-axis on plots
plot_order1 = [1:5,10,11,6:9,12:14];  
% Indices for Multiple-demand ROIs (not included in any of our main 
% analyses, but can plot results for these separately if you wish).
% the last three areas in the list are merged subregions of IPS (either all
% subregions together or two subregions at a time). 
plot_order2 = [15:20,22:24];

vismotor_names = ROI_names(plot_order1);
md_names = ROI_names(plot_order2);
plot_order_all = [plot_order1,plot_order2];
nROIs = length(plot_order_all);
vismotor_inds = find(ismember(plot_order_all,plot_order1));
md_inds = find(ismember(plot_order_all,plot_order2));

nVox2Use = 10000;
nPermIter=1000;

class_str = 'normEucDist';

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];

% all the different ways to train/test the decoder
condLabStrs = {'Trn/Test Predictable','Trn/Test Random','Trn Random/Test Predictable','Trn Predictable/Test Random'};
nConds = 2;

chance_val=0.5;

plotVisAcc = 1; % plot all four trn/test combinations for each area?
plotVisAccAvg=1; % average over both within-cond and both across-cond schemes?
plotMDAcc=1;
plotMDAccAvg=1;

diff_col=[0.5, 0.5, 0.5];
ms=10;  % marker size for significance dots
%% load results
nTrialsTotal = 400;

acc_allsubs = nan(nSubj,nROIs,nConds,2);    % last dim is trained on data from same or opp condition
d_allsubs = nan(nSubj,nROIs,nConds,2);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));

    save_dir = fullfile(curr_dir,'Decoding_results');

    % within-cond decoder
    fn2load = fullfile(save_dir,sprintf('TrnWithinCond_leavePairOut_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==length(ROI_names))
    acc_allsubs(ss,:,:,1) = mean(squeeze(allacc(plot_order_all,:,:)),3);
    d_allsubs(ss,:,:,1) = mean(squeeze(alld(plot_order_all,:,:)),3);
   
    % across-cond decoder
    fn2load = fullfile(save_dir,sprintf('TrnAcrossCond_leavePairOut_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==length(ROI_names))
    acc_allsubs(ss,:,:,2) = mean(squeeze(allacc(plot_order_all,:,:)),3);
    d_allsubs(ss,:,:,2) = mean(squeeze(alld(plot_order_all,:,:)),3);
  
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

% get some basic stats to use for the plots and tests below
vals = acc_allsubs;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));

%% make a bar plot of acc - visual areas
col = gray(6);
col= col(1:4,:);
if plotVisAcc
    
    mean_plot = [meanvals(vismotor_inds,:,1),meanvals(vismotor_inds,:,2)];
    sem_plot = [semvals(vismotor_inds,:,1),semvals(vismotor_inds,:,2)];
   
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,acclims,vismotor_names,condLabStrs,...
        'Accuracy','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
    
end

%%
col = gray(3);
col= col(1:2,:);

if plotVisAccAvg
    
    vals = mean(acc_allsubs(:,vismotor_inds,:,:),3);
    mean_plot = squeeze(mean(vals,1));
    sem_plot = squeeze(std(vals,[],1)./sqrt(nSubj));
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,acclims,vismotor_names,{'Within','Across'},...
        'Accuracy','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
    
end

%% make a bar plot of acc - md areas
col = gray(6);
col= col(1:4,:);
if plotMDAcc
    
    mean_plot = [meanvals(md_inds,:,1),meanvals(md_inds,:,2)];
    sem_plot = [semvals(md_inds,:,1),semvals(md_inds,:,2)];
   
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,acclims,md_names,condLabStrs,...
        'Accuracy','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
end

%%
col = gray(3);
col= col(1:2,:);

if plotMDAccAvg
    vals = mean(acc_allsubs(:,md_inds,:,:),3);
    mean_plot = squeeze(mean(vals,1));
    sem_plot = squeeze(std(vals,[],1)./sqrt(nSubj));
    plot_barsAndStars(mean_plot,sem_plot,[],...
        [],chance_val,acclims,md_names,{'Within','Across'},...
        'Accuracy','Spatial Memory Position',col)
    set(gcf,'Position',[800,800,1200,420])
end