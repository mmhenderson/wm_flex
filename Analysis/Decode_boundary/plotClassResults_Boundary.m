%% Plot accuracy of boundary orientation
% Decoding analysis itself performed in Classify_Boundary.m and saved as mat
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

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

% Indices into "ROI_names" corresponding to visual ROIs and motor ROIs
% reordering them a little for logical x-axis on plots
plot_order1 = [1:5,10,11,6:9,12:14];  
% Indices for Multiple-demand ROIs (not included in any of our main 
% analyses, but can plot results for these separately if you wish).
plot_order2 = [15:20]; 

vismotor_names = ROI_names(plot_order1);
md_names = ROI_names(plot_order2);
plot_order_all = [plot_order1, plot_order2];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
md_inds = find(ismember(plot_order_all,plot_order2));

nROIs = length(plot_order_all);

nVox2Use = 10000;
% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];
col = [125, 93, 175; 15, 127, 98]./255;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

chance_val=0.5;

plotVisMotorAcc = 1;    % make plots for retinotopic and motor ROIs?
plotMDAcc=0;    % make plots for MD ROIs?

%% load results

acc_allsubs = nan(nSubj,nROIs,nConds);
d_allsubs = nan(nSubj,nROIs,nConds);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyBoundary_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    acc_allsubs(ss,:,:) = mean(squeeze(allacc(plot_order_all,:,:)),3);
    d_allsubs(ss,:,:) = mean(squeeze(alld(plot_order_all,:,:)),3);
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

%% make a bar plot of acc - visual areas
if plotVisMotorAcc
   
    vals = squeeze(acc_allsubs(:,vismotor_inds,:));
    if nSubj>1
        meanvals =squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
    plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,vismotor_names,condLabStrs,'Accuracy','Classify Boundary Preview Orientation',col)
end


%% make a bar plot of acc - motor areas
if plotMDAcc
    
    vals = squeeze(acc_allsubs(:,md_inds,:));
    if nSubj>1
        meanvals = squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
    plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,md_names,condLabStrs,'Accuracy','Classify Boundary Preview Orientation',col)
end
