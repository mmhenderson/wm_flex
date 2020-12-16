%% plot signal related to contra/ipsi presses in digit localizer task
% this is the taks used to localize sensorimotor cortex areas.
% sanity check that these areas show the signal difference expected - this
% is a check for the response/finger labels.
%%
clear
close all;

sublist = [2:7];
nSubj=length(sublist);
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

condLabStrs = {'DigLocTask'};
nConds = length(condLabStrs);
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

nROIs = length(ROI_names);

vismotor_inds = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(vismotor_inds);
md_inds = [15:20];   % MD and motor ROIs
md_names = ROI_names(md_inds);

ylims = [-0.25, 0.25];
plotVisMotor=1;
plotVisMotorSS=1;
plotMD=0;

hemi_names={'LH','RH'};
resp_names={'Left Index', 'Right Index'};

nHemis=2;
nResp=2;

avg_sig = zeros(nSubj, nROIs, nHemis, nResp);

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('DigLocSignalByTrial_SepHemis_%s.mat',substr));
    load(fn2load);
    
    for vv = 1:nROIs
        for hh=1:nHemis
        
            %% pull out the data for main task

            if length(locSig)<vv || isempty(locSig(vv,hh).dat_avg) || size(locSig(vv,hh).dat_avg,2)<1
                fprintf('skipping %s area %s because not enough voxels\n',hemi_names{hh},ROI_names{vv})
                continue
            end

            for rr=1:nResp

                respLabs = locSig(vv,hh).ActDigit;
    
                digDat = locSig(vv,hh).dat_avg(respLabs==rr,:);
                
                avg_sig(ss,vv,hh,rr) = mean(digDat(:));
                
            end
        end
    end
end

%%
col = plasma(3);
col = col(1:2,:);
bw=0.50;
fs=14;
chance_val=0;
if plotVisMotor

    for hh=1:nHemis
        vals = squeeze(avg_sig(:,vismotor_inds,hh,:));
        if nSubj>1
            meanVals = squeeze(nanmean(vals,1));
            seVals = squeeze(nanstd(vals,[],1)./sqrt(nSubj));
        else
            meanVals = vals;
            seVals =[];
        end

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
        for cc=1:nResp
            b(cc).FaceColor = col(cc,:);
            b(cc).EdgeColor = col(cc,:);
            errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
                    'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
        end

        set(gca,'XTick', 1:numel(vismotor_inds))
        set(gca,'XTickLabel', vismotor_names,'XTickLabelRotation',90);
        ylabel('Mean Z-scored BOLD')
        set(gca,'YLim',ylims)
        set(gca,'XLim',[0,numel(vismotor_inds)+1])
        if chance_val~=0
            line([0,numel(vismotor_inds)+1],[chance_val,chance_val],'Color','k');
        end
        set(gca,'FontSize',fs);
        set(gcf,'Position',[800,800,1200,500]);
        % get locations of bars w offsets
        c=get(gcf,'Children');b=get(c(end),'Children');

        verspacerbig = range(ylims)/50;
        horspacer = abs(diff(bar_offset))/2;
    %     
        for vv=1:numel(vismotor_inds)
            % add individual subjects
            for ss=1:nSubj
                subvals = squeeze(avg_sig(ss,vismotor_inds(vv),hh,:));
                h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(5,:),'LineWidth',1.5);
                uistack(h,'bottom');
            end

        end
        b(end).BarWidth=bw;
        b(end-1).BarWidth=bw;
        leg=legend(lh,resp_names,'Location','EastOutside');
        set(gcf,'color','white')
        set(gcf, 'WindowStyle','normal','WindowState','normal')
        title(sprintf('%s: Mean BOLD during two conds of digit loc task',hemi_names{hh}));
        
    end
    

end

%%
if plotVisMotorSS

    motor_inds = 12:14;
    for hh=1:nHemis
        vals = squeeze(avg_sig(:,motor_inds,hh,:));
        if nSubj>1
            meanVals = squeeze(nanmean(vals,1));
            seVals = squeeze(nanstd(vals,[],1)./sqrt(nSubj));
        else
            meanVals = vals;
            seVals =[];
        end

        sub_colors = viridis(nSubj+1);
        set(groot,'DefaultLegendAutoUpdate','off');
        fh = figure();hold on;
        % first make the actual bar plot
        b = bar(gca,meanVals);
%         lh=[b(1),b(2)];

        % have to set this to "modal", otherwise it fails to get the XOffset
        % property.
        set(fh, 'WindowStyle','modal','WindowState','minimized')
        bar_offset = [b.XOffset];
        barPos = repmat((1:size(meanVals,1))', 1, length(bar_offset)) + repmat(bar_offset, size(meanVals,1), 1);
        for cc=1:nResp
            b(cc).FaceColor = col(cc,:);
            b(cc).EdgeColor = col(cc,:);
            b(cc).Visible='off';
%             errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
%                     'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
        end

        set(gca,'XTick', 1:numel(motor_inds))
        set(gca,'XTickLabel', vismotor_names(motor_inds),'XTickLabelRotation',90);
        ylabel('Mean Z-scored BOLD')
        set(gca,'YLim',ylims)
        set(gca,'XLim',[0,numel(motor_inds)+1])
        if chance_val~=0
            line([0,numel(motor_inds)+1],[chance_val,chance_val],'Color','k');
        end
        set(gca,'FontSize',fs);
        set(gcf,'Position',[800,800,1200,500]);
        % get locations of bars w offsets
        c=get(gcf,'Children');b=get(c(end),'Children');

        verspacerbig = range(ylims)/50;
        horspacer = abs(diff(bar_offset))/2;
    %     
        for vv=1:numel(motor_inds)
            % add individual subjects
            lh=[];
            for ss=1:nSubj
                subvals = squeeze(avg_sig(ss,motor_inds(vv),hh,:));
                h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(ss,:),'LineWidth',1.5);
                lh=[lh,h];
                ll{ss} = sprintf('S%02d',sublist(ss));
                uistack(h,'bottom');
            end

        end
        b(end).BarWidth=bw;
        b(end-1).BarWidth=bw;
        leg=legend(lh,ll,'Location','EastOutside');
        set(gcf,'color','white')
        set(gcf, 'WindowStyle','normal','WindowState','normal')
        title(sprintf('%s: Mean BOLD during two conds of digit loc task',hemi_names{hh}));
        
    end
    

end

%%
if plotMD
    
    for hh=1:nHemis
        vals = squeeze(avg_sig(:,md_inds,hh,:));
        if nSubj>1
            meanVals = squeeze(nanmean(vals,1));
            seVals = squeeze(nanstd(vals,[],1)./sqrt(nSubj));
        else
            meanVals = vals;
            seVals =[];
        end

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
        for cc=1:nResp
            b(cc).FaceColor = col(cc,:);
            b(cc).EdgeColor = col(cc,:);
            errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
                    'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
        end

        set(gca,'XTick', 1:numel(md_inds))
        set(gca,'XTickLabel', md_names,'XTickLabelRotation',90);
        ylabel('Mean Z-scored BOLD')
        set(gca,'YLim',ylims)
        set(gca,'XLim',[0,numel(md_inds)+1])
        if chance_val~=0
            line([0,numel(md_inds)+1],[chance_val,chance_val],'Color','k');
        end
        set(gca,'FontSize',fs);
        set(gcf,'Position',[800,800,1200,500]);
        % get locations of bars w offsets
        c=get(gcf,'Children');b=get(c(end),'Children');

        verspacerbig = range(ylims)/50;
        horspacer = abs(diff(bar_offset))/2;
    %     
        for vv=1:numel(md_inds)
            % add individual subjects
            for ss=1:nSubj
                subvals = squeeze(avg_sig(ss,md_inds(vv),hh,:));
                h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(5,:),'LineWidth',1.5);
                uistack(h,'bottom');
            end

        end
        b(end).BarWidth=bw;
        b(end-1).BarWidth=bw;
        leg=legend(lh,resp_names,'Location','EastOutside');
        set(gcf,'color','white')
        set(gcf, 'WindowStyle','normal','WindowState','normal')
        title(sprintf('%s: Mean BOLD during two conds of digit loc task',hemi_names{hh}));
        
    end
    

end