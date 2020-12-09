% script to plot the result of decoding analyses for oriSpin. 

clear
close all;

sublist = [2,3,5,6];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;

% names of the ROIs 
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S postcentral','G postcentral','S central','G precentral',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
nROIs = length(ROI_names);

% plotVisD = 1;
plotVisualAcc = 1;
% plotMotorD = 1;
plotMotorMDAcc = 1;

plot_order1 = [1:5,10,11,6:9];  % visual ROIs
visual_names = ROI_names(plot_order1);
plot_order2 = [12:21];   % MD and motor ROIs
motor_names = ROI_names(plot_order2);
if ~plotVisualAcc
    plot_order1 = [];
end
plot_order_all = [plot_order1, plot_order2];
vis_inds = find(ismember(plot_order_all,plot_order1));
motor_inds = find(ismember(plot_order_all,plot_order2));

nVOIs = length(plot_order_all);

nVox2Use = 10000;
% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4,0.75];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);
% cc=1;

% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];
chance_val=0.5;

%% load results
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
lw =1;


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);


acc_allsubs = nan(nSubj,nVOIs,nConds,nTRs_out);
d_allsubs = nan(nSubj,nVOIs,nConds,nTRs_out);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyDiskSides_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    
    acc_allsubs(ss,:,:,:) = mean(squeeze(allacc(plot_order_all,:,:,:)),3);
    d_allsubs(ss,:,:,:) = mean(squeeze(alld(plot_order_all,:,:,:)),3);
    
end

assert(~any(isnan(acc_allsubs(:))))
assert(~any(isnan(d_allsubs(:))))



%% Now we make a figure for all subjects
cc=1;
for vi = 1:nVOIs
   
    if plotVisualAcc && ismember(vi, vis_inds)
        nplots_vis = ceil(sqrt(numel(vis_inds)));
        
        figure(1);hold all;
        subplot(nplots_vis,ceil(numel(vis_inds)/nplots_vis),vi);hold all;
        for cc = 1:nConds
            vals = acc_allsubs(:,vi,cc,:);
            if nSubj==1
                meanvals =squeeze(vals);
                semvals = [];
            else
                meanvals = squeeze(mean(vals,1));
                semvals = squeeze(std(vals,[],1))./sqrt(nSubj);
            end
        
            errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
        end
        
        set(gca, 'FontSize', 12, 'XLim',[0 max(tax)])
        ylim(acclims);
        plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
        end
        
        if vi==1
            xlabel('Time(s)');
            ylabel('Accuracyy');
        end

        if vi==numel(vis_inds)
            legend(condLabStrs);
        end
        
        if contains(visual_names{vi}, ' ')
            % break it into two strings
            spaceind = find(visual_names{vi}==' ');
            title(sprintf('%s\n%s', visual_names{vi}(1:spaceind-1), visual_names{vi}(spaceind+1:end)));
        else
            title(sprintf('%s', visual_names{vi}));
        end
       
      
    end

    if plotMotorMDAcc && ismember(vi,motor_inds)
        nplots_motor = ceil(sqrt(numel(motor_inds)));
        
        figure(2);hold all;

        subplot(nplots_motor,ceil(numel(motor_inds)/nplots_motor),vi-numel(vis_inds));hold all;

        for cc = 1:nConds
            vals = acc_allsubs(:,vi,cc,:);
            if nSubj==1
                meanvals =squeeze(vals);
                semvals = [];
            else
                meanvals = squeeze(mean(vals,1));
                semvals = squeeze(std(vals,[],1))./sqrt(nSubj);
            end

            errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
        end
        
        set(gca, 'FontSize', 12, 'XLim',[0 max(tax)])
        ylim(acclims);
        plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
        end
        if vi-numel(vis_inds)==1
            xlabel('Time(s)');
            ylabel('Accuracy');
        end

        if vi==numel(plot_order_all)
            legend(condLabStrs);
        end
        
        if contains(motor_names{vi-numel(vis_inds)}, ' ')
            % break it into two strings
            spaceind = find(motor_names{vi-numel(vis_inds)}==' ');
            title(sprintf('%s\n%s', motor_names{vi-numel(vis_inds)}(1:spaceind-1), motor_names{vi-numel(vis_inds)}(spaceind+1:end)));
        else
            title(sprintf('%s', motor_names{vi-numel(vis_inds)}));
        end
       
    end
    
end
%%
if plotVisualAcc
    figure(1);hold all;
%     suptitle('Fidelity Over Time: Train Localizer');
     set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
end
if plotMotorMDAcc
    figure(2);hold all;
%     suptitle('Fidelity Over Time: Train Localizer');
     set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
end