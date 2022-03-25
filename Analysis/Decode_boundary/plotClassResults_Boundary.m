%% Plot accuracy of boundary orientation
% Decoding analysis itself performed in Classify_Boundary.m and saved as mat
% file. This script loads that file, does all stats and plotting. 
% note this analysis is not in our paper.
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

condLabStrs = {'Informative','Uninformative'};
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
    fn2load = fullfile(save_dir,sprintf('ClassifyBoundary360_%s_%dvox_%s.mat',class_str,nVox2Use,substr));    
%     fn2load = fullfile(save_dir,sprintf('ClassifyBoundary_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    acc_allsubs(ss,:,:) = mean(squeeze(allacc(plot_order_all,:,:)),3);
    d_allsubs(ss,:,:) = mean(squeeze(alld(plot_order_all,:,:)),3);
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))

%% make a bar plot of decoding acc, with single subjects overlaid
bw=0.50;
fs=14;

if plotVisMotorAcc
    
    vals = squeeze(acc_allsubs(:,vismotor_inds,:));
    if nSubj>1
        meanvals =squeeze(mean(vals,1));
        semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    else
        meanvals = vals;
        semvals =[];
    end
   
    meanVals=meanvals;
    seVals=semvals;
    
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
%         % add significance of individual areas/conditions
%         for cc=1:nConds
%             for aa=1:numel(alpha_vals)
%                 if p_sr(vismotor_inds(vv),cc)<alpha_vals(aa)
%                     % smaller dots get over-drawn with larger dots
%                     plot(vv+bar_offset(cc), meanVals(vv,cc)+seVals(vv,cc)+verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa))
%                 end
%             end
%         end
%         % add significance of condition differences
%         for aa=1:numel(alpha_vals)
%             if p_diff(vismotor_inds(vv))<alpha_vals(aa)
%                 [mx,maxind] = max(meanVals(vv,:));
%                 % smaller dots get over-drawn with larger dots
%                 plot(vv+bar_offset, repmat(meanVals(vv,maxind)+seVals(vv,maxind)+2*verspacerbig,2,1),'-','Color','k','LineWidth',1)
%                 plot(vv, meanVals(vv,maxind)+seVals(vv,maxind)+3*verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa));
%                 
%             end
%             if vv==1
%                 lh=[lh,plot(-1, meanVals(vv,1)+seVals(vv,1)+3*verspacerbig,'.','Color','k','MarkerSize',alpha_ms(aa))];
%             end
%         end
    end
    b(end).BarWidth=bw;
    b(end-1).BarWidth=bw;
    legend(lh, {'Informative', 'Uninformative'})
%     leg=legend(lh,{'Predictable','Random','p<0.05','0<0.01','p<0.001'},'Location','EastOutside');

    set(gcf,'color','white')
    set(gcf, 'WindowStyle','normal','WindowState','normal')
    title('Classify orientation of "preview" disk (use 0-360 space)')
%     title('Classify orientation of "preview" disk (use 0-180 space)')
%     saveas(gcf,fullfile(figpath,'TrainSWM_TestWM_allareas.pdf'),'pdf');
end
% 
% %% make a bar plot of acc - visual areas
% fs=14;
% if plotVisMotorAcc
%    
%     vals = squeeze(acc_allsubs(:,vismotor_inds,:));
%     if nSubj>1
%         meanvals =squeeze(mean(vals,1));
%         semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
%     else
%         meanvals = vals;
%         semvals =[];
%     end
%     plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,vismotor_names,condLabStrs,...
%         'Accuracy','Classify orientation of "preview" disk (use 0-360 space)',col)
%     set(gca,'FontSize',fs);
%     set(gcf,'Position',[800,800,1200,500]);
% end
% 
% 
% %% make a bar plot of acc - motor areas
% if plotMDAcc
%     
%     vals = squeeze(acc_allsubs(:,md_inds,:));
%     if nSubj>1
%         meanvals = squeeze(mean(vals,1));
%         semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
%     else
%         meanvals = vals;
%         semvals =[];
%     end
%     plot_barsAndStars(meanvals,semvals,[],[],chance_val,acclims,md_names,condLabStrs,'Accuracy',...
%         'Classify orientation of "preview" disk (use 0-360 space)',col)
% end
