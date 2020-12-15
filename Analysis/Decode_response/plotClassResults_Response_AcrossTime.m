% script to plot the result of decoding analyses for oriSpin. 

clear
close all;

sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
     'S1','M1','Premotor',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1-S1 all'};

plot_order = [12:14];  % visual ROIs
% plot_order=[4,5];
plot_names = ROI_names(plot_order);

nROIs = length(plot_order);

nVox2Use = 10000;
nPermIter=1000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4,1];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);
% cc=1;

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];
chance_val=0.5;

sig_heights = [0.91,0.93,0.97];
diff_col=[0.5, 0.5, 0.5];

fs=12;  % font size for all plots
ms=10;  % marker size for significance dots
%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
lw =1;


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);


acc_allsubs = nan(nSubj,nROIs,nConds,nTRs_out, nTRs_out);
d_allsubs = nan(nSubj,nROIs,nConds,nTRs_out, nTRs_out);
% accrand_allsubs = nan(nSubj, nVOIs, nConds, nTRs_out, nPermIter);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
%     fn2load = fullfile(save_dir,sprintf('TrnWithinCond_AcrossTime_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_AcrossTime_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names))
    acc_allsubs(ss,:,:,:,:) = squeeze(allacc(plot_order,:,:,:));
    d_allsubs(ss,:,:,:,:) = squeeze(alld(plot_order,:,:,:));
    
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))


%% Now make a plot for each ROI of interest

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
    % 7 is the first condition diagonals
    % 6 is the first condition matrix
    % 4 is the second condition diagonals
    % 2 is the second condition matrix
%     c(7).Position(1) = c(6).Position(1);
%     c(7).Position(3) = c(6).Position(3);
    c(7).OuterPosition(1) = c(6).OuterPosition(1);
    c(7).OuterPosition(3) = c(6).OuterPosition(3);
%     c(4).Position(1) = c(2).Position(1);  
%     c(4).Position(3) = c(2).Position(3);
    c(4).OuterPosition(1) = c(2).OuterPosition(1);
    c(4).OuterPosition(3) = c(2).OuterPosition(3);
end