%% Plot average signal in each ROI over time. 
% similar to deconvolution but instead of deconvolving, just stack epoched
% trials and take average. 

%%
clear
close all;
sublist = [2:7];
nSubj = length(sublist);
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
figpath = fullfile(exp_path,'figs');
addpath(fullfile(exp_path,'Analysis','stats_code'));

% names of the ROIs - note the motor areas have been changed to nicer names
% here for plotting.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc'...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
nROIs = length(plot_order_all);

lw=1;

% tr2baseline=0;

condLabStrs = {'Predictable','Random'};

% length of each run type
% nTRs = 583-16;
% trDur = 0.80;

% deconvolution param
% nTRs_out = 30; % # of TRs post-stimulus to model in the decon analysis
% define timing axis
% tax = trDur*(0:nTRs_out-1); 

nTrialsPerRun = 20;
nConds = 2;

allsub_bold_mean = nan(nSubj, length(plot_order_all), nConds);
allsub_bold_sem = nan(nSubj, length(plot_order_all), nConds);

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);

% % events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

ylims=[-0.25, 0.25];
% sig_heights = [0.40,0.44,0.48];
avg_TRs_range = [10,16];

col_conds = [125, 93, 175; 15, 127, 98]./255;

diff_col=[0.5, 0.5, 0.5];
%%
for ss=1:length(sublist)

    
    substr = sprintf('S%02d',sublist(ss));
   
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    
    %
    for vi = 1:nROIs % for all visual areas I want to look at
        vv = plot_order_all(vi);
        
        % trials by TRs by voxels
        dat_avg = mean(mainSig(vv).dat_by_TR(:,avg_TRs_range(1)+1:avg_TRs_range(2)+1,:),2);
        assert(all(all(dat_avg==mainSig(vv).dat_avg)));
        % baseline to zero at first TR
%         dat_by_TR = dat_by_TR - repmat(dat_by_TR(:,1,:),1,size(dat_by_TR,2),1);
        
        condLabs = mainSig(vv).condLabs;
        nVox = size(dat_avg,2);
        
        % take an average over voxels and save it for later
        for cc = 1:nConds
           mean_over_trials= squeeze(mean(dat_avg(condLabs==cc,:),1));
           allsub_bold_mean(ss,vi,cc) = mean(mean_over_trials,2);
           allsub_bold_sem(ss,vi,cc) = std(mean_over_trials,[],2)./sqrt(size(mean_over_trials,2));
        end
       
        if nVox==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end

        fprintf('processing %s area %s, %d voxels\n', substr, ROI_names{vv}, nVox);
        
      
 
    end

end


%% compare between conditions
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 456556;
rng(rndseed,'twister')
nPermIter=1000;
vals = allsub_bold_mean;
real_sr_stat = nan(nROIs,1);
rand_sr_stat = nan(nROIs, nPermIter);

% p_diff_sr=nan(nROIs,nTRs_out);
for vv=1:nROIs
   
        realvals = squeeze(vals(:,vv,:));
        % what is the sign-rank statistic for the real data?
        real_sr_stat(vv) = signrank_MMH(realvals(:,1),realvals(:,2));
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
p_cond_diff = 2*min(cat(2,mean(repmat(real_sr_stat,1,nPermIter)>=rand_sr_stat,2), ...
    mean(repmat(real_sr_stat,1,nPermIter)<=rand_sr_stat,2)),[],2);

% array2table(squeeze(p_cond_diff(vismotor_inds,1:15)),...
%     'RowNames',vismotor_names,'VariableNames',strseq('CondDiff_TR',1:15))
% array2table(squeeze(p_cond_diff(vismotor_inds,16:30)),...
%     'RowNames',vismotor_names,'VariableNames',strseq('CondDiff_TR',16:30))

%%
fs=14;
bw=0.50;
% if plotVisMotor
    
    % get some basic stats to use for the plots and tests below
    vals = allsub_bold_mean;
    meanvals = squeeze(mean(vals,1));
    semvals = squeeze(std(vals,[],1)./sqrt(nSubj));

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
        b(cc).FaceColor = col_conds(cc,:);
        b(cc).EdgeColor = col_conds(cc,:);
        errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    end

    set(gca,'XTick', 1:numel(vismotor_inds))
    set(gca,'XTickLabel', vismotor_names,'XTickLabelRotation',90);
    ylabel('BOLD (z-score)')
    set(gca,'YLim',ylims)
    set(gca,'XLim',[0,numel(vismotor_inds)+1])
%     if chance_val~=0
%         line([0,numel(vismotor_inds)+1],[chance_val,chance_val],'Color','k');
%     end
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
            subvals = squeeze(vals(ss,vismotor_inds(vv),:));
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
        % add significance of condition differences
        for aa=1:numel(alpha_vals)
            if p_cond_diff(vismotor_inds(vv))<alpha_vals(aa)
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
   
    saveas(gcf,fullfile(figpath,'Univariate_sepconds_avgtime_allareas.pdf'),'pdf');

    
% end

