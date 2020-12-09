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
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

nROIs = length(ROI_names);

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vimotor_names = ROI_names(plot_order1);
plot_order2 = [15:20];   % MD and motor ROIs
md_names = ROI_names(plot_order2);
plot_order_all = [plot_order1, plot_order2];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
md_inds = find(ismember(plot_order_all,plot_order2));

nVOIs = length(plot_order_all);

% nVox2Use = 100;
nVox2Use=10000;
% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

chance_val=0.5;

plotVisMotorD = 0;
plotVisMotorAcc = 1;
plotMD_D = 0;
plotMD_Acc = 1;
%% load results

acc_allsubs = nan(nSubj,nVOIs,nConds);
d_allsubs = nan(nSubj,nVOIs,nConds);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyDiskSides_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
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
    plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,vimotor_names,condLabStrs,'Accuracy','Classify Disk Sides (average of four 2-way)',col)
end
%% make a bar plot of dprimean(me - visual areas
if plotVisMotorD
    
    vals = squeeze(d_allsubs(:,vismotor_inds,:));
    if nSubj>1
        meanvals =squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
    plot_barsAndStars(meanvals,semvals,[],[],[],dprimelims,vimotor_names,condLabStrs,'Dprime','Classify Disk Sides (average of four 2-way)',col)
end

%% make a bar plot of acc - motor areas
if plotMD_Acc
    
    vals = squeeze(acc_allsubs(:,md_inds,:));
    if nSubj>1
        meanvals =squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
    plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,md_names,condLabStrs,'Accuracy','Classify Disk Sides (average of four 2-way)',col)
end
%% make a bar plot of dprime - motor areas
if plotMD_D
   
    vals = squeeze(d_allsubs(:,md_inds,:));
    if nSubj>1
        meanvals =squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
    plot_barsAndStars(meanvals,semvals,[],[],[],dprimelims,md_names,condLabStrs,'Dprime','Classify Disk Sides (average of four 2-way)',col)
end