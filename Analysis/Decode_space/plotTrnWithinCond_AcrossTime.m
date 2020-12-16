%% Plot accuracy of spatial decoder
% trained within each condition of main task. Trained and tested across
% timepoints - so get a full [nTR x nTR] temporal generalization matrix.
% Decoding analysis itself performed in TrnWithinCond_AcrossTime_leavePairOut.m
% and saved as mat file. 
% This script loads that file, does all stats and plotting. 
%%
clear
close all;

sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;

% names of the ROIs 
% the last three areas in the list are merged subregions of IPS (either all
% subregions together or two subregions at a time). 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
     'S1','M1','Premotor',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1-S1 all',...
    'IPS0-3','IPS0-1','IPS2-3'};

% Indices into "ROI_names" corresponding to visual ROIs and motor ROIs
plot_order = [1:5,10,11,6:9,22:24];  % visual ROIs
plot_names = ROI_names(plot_order);

nROIs = length(plot_order);

nVox2Use = 10000;
nPermIter=1000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4,1];
dprimelims = [-0.2, 1.4];
col = [125, 93, 175; 15, 127, 98]./255;

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];
chance_val=0.5;

fs=12;  % font size for all plots
%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1); % time axis, how many seconds is each TR
lw =1;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

acc_allsubs = nan(nSubj,nROIs,nConds,nTRs_out, nTRs_out);
d_allsubs = nan(nSubj,nROIs,nConds,nTRs_out, nTRs_out);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('TrnWithinCond_AcrossTime_leavePairOut_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names))
    % averaging decoding performance over the four decoding schemes (0 vs
    % 180, 45 versus 225, etc)
    acc_allsubs(ss,:,:,:,:) = mean(squeeze(allacc(plot_order,:,:,:,:)),3);
    d_allsubs(ss,:,:,:,:) = mean(squeeze(alld(plot_order,:,:,:,:)),3);    
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))


%% Now make a plot for each ROI

my_clim=[0.33, 0.865];
minval=1;
maxval=0;
step=5;
for vi = 1:numel(plot_order)
    
    figure();hold all;
    colormap(plasma);

    pis = [1,3;2,4];
    for cc=1:nConds

        % first plot the diagonal only
        subplot(2,2,pis(cc,1));hold all;
        vals = squeeze(acc_allsubs(:,vi,cc,:,:));
        if nSubj==1
            meanvals =squeeze(vals)';
            semvals = [];
        else
            meanvals = squeeze(mean(vals,1));
            semvals = squeeze(std(vals,[],1))./sqrt(nSubj);
        end
        plot(tax,diag(meanvals),'-','Color',col(cc,:),'LineWidth',lw);
        bandedError_MMH(tax, diag(meanvals)',diag(semvals)', col(cc,:), 0.5);
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)]);
        ylim(acclims);
        plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
        end       
        xlabel('Time(s)');
        ylabel('Accuracy');
        
        title(sprintf('%s', condLabStrs{cc}));
        
        % now plot the whole conf matrix
        subplot(2,2,pis(cc,2));hold all;
       
        minval=min([min(meanvals(:)), minval]);
        maxval=max([max(meanvals(:)), maxval]);
        assert(minval>my_clim(1) & maxval<my_clim(2));
        
        imagesc(meanvals,[my_clim]);
        colorbar()
        axis square equal 
        set(gca,'XTick',[1:step:size(meanvals,1)],'XTickLabel',tax(1:step:end),'YTick',[1:step:size(meanvals,1)],'YTickLabel',tax(1:step:end))
        
        % adding lines to this plot to indicate the start and end of trial
        % events (note these won't line up with the imagesc grid, because
        % events don't fall on exact TR boundaries.)
        for ee=1:numel(evts2plot)
            % convert from seconds into TRs here, so we can add to the plot
            % in its original coordinate system (ticks are TRs)
            % since TR 1 on this plot is 0 seconds, need to add a 1 as
            % well.
            loc = evts2plot(ee)/trDur+1;
            plot([loc,loc],get(gca,'YLim'),'-','Color','k')
            plot(get(gca,'XLim'),[loc,loc],'-','Color','k')
            
        end      
        set(gca, 'FontSize', fs, 'XLim',[0.5, size(meanvals,1)+0.5],'YLim',[0.5, size(meanvals,1)+0.5])        
        xlabel('Testing timept (s)');
        ylabel('Training timept (s)');

    end
    suptitle(sprintf('%s\nTrain within condition',plot_names{vi}));
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800])
    c=get(gcf,'Children');  
    % trying to make subplots line up in width...
    % 7 is the first condition diagonals
    % 6 is the first condition matrix
    % 4 is the second condition diagonals
    % 2 is the second condition matrix
    c(7).OuterPosition(1) = c(6).OuterPosition(1);
    c(7).OuterPosition(3) = c(6).OuterPosition(3);
    c(4).OuterPosition(1) = c(2).OuterPosition(1);
    c(4).OuterPosition(3) = c(2).OuterPosition(3);
end