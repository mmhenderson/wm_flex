% plot reconstructions 

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
addpath(fullfile(exp_path,'Analysis','stats_code'))

% names of the ROIs - note the motor areas have been changed to nicer names
% here for plotting.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

plot_order1 = [1:5,10,11,6:9,12:14];  % visual ROIs
vismotor_names = ROI_names(plot_order1);
plot_order2 = [15:20];   % MD and motor ROIs
md_names = ROI_names(plot_order2);
plot_order_all = [plot_order1, plot_order2];
vismotor_inds = find(ismember(plot_order_all,plot_order1));
md_inds = find(ismember(plot_order_all,plot_order2));

nROIs = length(plot_order_all);

plotVisMotorRecons = 1;
plotVisMotorReconsAvg=1;
plotMDRecons = 0;
plotMDReconsAvg=0;

plotVisMotorFids = 1;
plotVisualFidSideMatch = 0;
plotMDFids = 1;
plotMDFidSideMatch=0;

plotVisBiasSepDirs = 0;
plotVisBiasAvgDirs = 0;
plotVisualBiasBars=0;
plotMotorBiasBars=0;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

nVOIs = length(plot_order_all);

biaslims=[-180,180];
ylims = [-0.4,0.4];
fidlims = [-0.1, 0.25];
% col = plasma(5);
% col = col(2:2:end-1
col = [125, 93, 175; 15, 127, 98]./255;

fs=20;

alpha_vals=[0.05, 0.01, 0.001];
alpha_ms = [8,16,24];
alpha = alpha_vals(1);

col_resp = plasma(3);
col_resp = col_resp(1:2,:);

% when looking at category repulsion effects, how many distance bins to use?
nBins=1;
dist_bin_edges = linspace(0,90,nBins+1);
for bb=1:nBins
    binLabs{bb} = sprintf('%d-%d deg from bound', dist_bin_edges(bb),dist_bin_edges(bb+1));
end
nDirs =  2; % CW/CCW
dirLabs = {'Targ CW (-) of bound','Targ CCW (+) of bound'};

sameSideLabs = {'Finger on same side as target','Finger on opp side of target'};

%%
for ss=1:nSubj
    
    substr = sprintf('S%02d',sublist(ss));    
    fn = fullfile(exp_path,'Analysis','IEM','IEM_results',sprintf('TrnSWMLoc_TestWM_%s.mat',substr));
    load(fn);
    assert(numel(allchanresp)==numel(ROI_names));
    if ss==1        
        %preallocate some arrays       
        avg_recs = nan(nSubj, nVOIs, nConds, length(xx));
        recs_binned_preview_bound = nan(nSubj, nVOIs, nConds, nBins, nDirs, length(xx));
        
        fidelity = nan(nSubj, nVOIs, nConds);
        fidelity_side_match = nan(nSubj, nVOIs, nConds, length(sameSideLabs));
        bias = nan(nSubj, nVOIs, nConds, nBins);
    end
    
    condlabs = allchanresp(1).condLabs;

    dist_to_bound = allchanresp(1).dist_to_real_bound;
    which_dir = allchanresp(1).dir_to_real_bound;
    % for random condition, using the random stimulus that was shown to
    % compute bias from the boundary (otherwise this is meaningless)
    dist_to_bound(condlabs==2) = allchanresp(1).dist_to_rand_bound(condlabs==2);
    which_dir(condlabs==2) = allchanresp(1).dir_to_rand_bound(condlabs==2);
     
    % calculate what the response should be for the RANDOM TRIALS
    % at this point (what is the response associated with the random
    % preview disk)? this is irrelevant to the task, but may be 
    % automatically represented.
    targ_pos = allchanresp(1).targPos;
    preview_bound_pos = allchanresp(1).randBoundPos;
    preview_bound_pos(condlabs==1) = allchanresp(1).boundPos(condlabs==1);
    resp_assoc_preview = zeros(size(targ_pos));
    over180 = preview_bound_pos>180;
    resp_assoc_preview(over180 &...
        targ_pos<preview_bound_pos &...
        targ_pos>mod(preview_bound_pos+180,360)) = 2;
    resp_assoc_preview(over180 & ...
        (targ_pos>preview_bound_pos |...
        targ_pos<mod(preview_bound_pos+180,360))) = 1;
    resp_assoc_preview(~over180 &...
        targ_pos>preview_bound_pos &...
        targ_pos<mod(preview_bound_pos+180,360)) = 1;
    resp_assoc_preview(~over180 &...
        (targ_pos<preview_bound_pos |...
        targ_pos>mod(preview_bound_pos+180,360))) = 2;
    
    % is the target on the same half of the screen as the prepared
    % response? e.g. preparing a left finger response while maintaining a 
    % representation of item on left of screen?
    % resp 1=L, 2=R. 
    same_side = 2*ones(size(resp_assoc_preview));
    same_side((targ_pos<90 | targ_pos>270) & resp_assoc_preview==2) = 1; % R
    same_side((targ_pos>90 & targ_pos<270) & resp_assoc_preview==1) = 1; % L
    
    which_bin = zeros(size(dist_to_bound));
    for bb=1:nBins
        inds2use = (dist_to_bound<dist_bin_edges(bb+1) &...
                   dist_to_bound>dist_bin_edges(bb));
        which_bin(inds2use) = bb;
    end
    
    %% load recons from each area
    for vv=1:nVOIs
      
       if plot_order_all(vv)>length(allchanresp) || isempty(allchanresp(plot_order_all(vv)).chan_resp_shift) 
           fprintf('skipping %s for S%s because no voxels\n', ROI_names{plot_order_all(vv)}, substr);
           continue
       end
       
       
       for cc = 1:nConds
           
           % take out trials from one condition at a time
           theserecs = allchanresp(plot_order_all(vv)).chan_resp_shift(condlabs==cc,:);
           
           avg_recs(ss,vv,cc,:) = mean(theserecs,1); 

           
           for bb=1:nBins
               
               for ww=1:nDirs
               
                   % ww=2 for CCW, ww=1 for CW
                   inds2use = which_bin(condlabs==cc)==bb & which_dir(condlabs==cc)==ww;
                   recs_this_bin = theserecs(inds2use,:);
                   recs_binned_preview_bound(ss,vv,cc,bb,ww,:) = mean(recs_this_bin,1);
                
               end
               
               recs1 = squeeze(recs_binned_preview_bound(ss,vv,cc,bb,1,:));
               recs2 = squeeze(recs_binned_preview_bound(ss,vv,cc,bb,2,:));
               rec = mean(cat(2,flipud(recs1), recs2),2);
               [~,peakind] = max(rec);
               bias(ss,vv,cc,bb) = xx(peakind)-180;
           end
           
           % get the fidelity
           angs = abs((1:360)-180);
           cos_vals = cosd(angs);
           fidelity(ss,vv,cc) = mean(cos_vals.*mean(theserecs,1));

           % compute fidelity separately for trials where target was on
           % same side as response (e.g. left finger/left side of screen)
           for mm=1:2
               theserecs = allchanresp(plot_order_all(vv)).chan_resp_shift(condlabs==cc & same_side==mm,:);           
               fidelity_side_match(ss,vv,cc,mm) = mean(cos_vals.*mean(theserecs,1));
           end
           
           if ss>1 && fidelity(ss-1,vv,cc)==fidelity(ss,vv,cc)
               error('you have a missing value and something bad is happening')
           end
           
       end
    end
  
end

%% make plots: recons in visual ROIs

lh = [];
ylims_big =[-0.4, 0.6];
if plotVisMotorRecons
    figure;hold all;
    
    % plot each subject and ROI as a separate subplot.
    for ss = 1:nSubj
        for vv = 1:numel(vismotor_inds)
            
            subplot(nSubj+1,numel(vismotor_inds),(ss-1)*numel(vismotor_inds)+vv);hold all;
            for cc = 1:nConds
                plot(xx,squeeze(avg_recs(ss,vismotor_inds(vv),cc,:)),'Color',col(cc,:),'LineWidth',1);
                
            end
            set(gca, 'FontSize', 12)
            set(gca,'XTick',[]);
            set(gca,'YLim',ylims_big);
%             set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims)

            plot([shift_to shift_to], ylims_big, 'k', 'LineWidth', 1)
            if vv ==1
                ylabel(sprintf('S%02d',sublist(ss)));
            end

            if ss==1
                if contains(vismotor_names{vv}, ' ')
                    % break it into two strings 
                    spaceind = find(vismotor_names{vv}==' ');
                    title(sprintf('%s\n%s', vismotor_names{vv}(1:spaceind-1), vismotor_names{vv}(spaceind+1:end)));
                else               
                    title(sprintf('%s', vismotor_names{vv}));
                end
            end
            
        end
    end
    % now plot the average across subjects.
    for vv = 1:numel(vismotor_inds)
            
        subplot(nSubj+1,numel(vismotor_inds),(nSubj)*numel(vismotor_inds)+vv);hold all;
      
        for cc = 1:nConds
            recs = squeeze(avg_recs(:,vismotor_inds(vv),cc,:));
            if nSubj>1
                meanvals = nanmean(recs,1);
            else
                meanvals = recs;
            end
            h= errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);
            if vv==numel(vismotor_inds)
                lh = [lh,h];
            end
        end
        
        set(gca, 'FontSize', 12)
        if vv==1
            set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims)
        else
            set(gca,'XTick',[]);
            set(gca,'YLim',ylims,'YTick',[ylims(1),0,ylims(2)]);
        end
           
        plot([shift_to shift_to], ylims, 'k', 'LineWidth', 1)
        if vv ==1
            ylabel('Average');           
            ylabel('Avg');
            xlabel('Orientation Channel')    
        end
        
        set(gcf,'Color','w')
    end
    legend(lh,condLabStrs,'Location','SouthEast');
    suptitle('Train Loc, Test WM')
    set(gcf,'Position',[200,200,1400,800]);
end

%% subject average

if plotVisMotorReconsAvg
    lh = [];
    npx = ceil(sqrt(numel(vismotor_inds)));
    npy = ceil(numel(vismotor_inds)/npx);
    figure;hold all;
   
    % now plot the average across subjects.
    for vv = 1:numel(vismotor_inds)
            
        subplot(npx,npy,vv);hold all;
      
        for cc = 1:nConds
            recs = squeeze(avg_recs(:,vismotor_inds(vv),cc,:));
            if nSubj>1
                meanvals = nanmean(recs,1);
            else
                meanvals = recs;
            end
            h= errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);
            if vv==numel(vismotor_inds)
                lh = [lh,h];
            end
        end
        
        set(gca, 'FontSize', fs)
        if vv==numel(vismotor_inds)-1
            set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims,'YTick',[ylims(1),0,ylims(2)])
            ylabel('a.u.');
            xlabel('Spatial Position (deg)')    
        else
            set(gca,'XLim', [0, 360],'XTick',[]);
            set(gca,'YLim',ylims,'YTick',[ylims(1),0,ylims(2)]);
        end
        title(vismotor_names{vv})
        plot([shift_to shift_to], ylims, 'k', 'LineWidth', 1)
        
        set(gcf,'Color','w')
    end
    legend(lh,condLabStrs,'Location','SouthEast');
    suptitle('Train Loc, Test WM')
    set(gcf,'Position',[200,200,1400,800]);
end

%% make plots: recons in motor/MD ROIs
lh = [];
if plotMDRecons
    figure;hold all;
    
    % plot each subject and ROI as a separate subplot.
    for ss = 1:nSubj
        for vv = 1:numel(md_inds)
            
            subplot(nSubj+1,numel(md_inds),(ss-1)*numel(md_inds)+vv);hold all;
            for cc = 1:nConds
                plot(xx,squeeze(avg_recs(ss,md_inds(vv),cc,:)),'Color',col(cc,:),'LineWidth',1);
                
            end
            set(gca, 'FontSize', 12)
            set(gca,'XTick',[]);
            set(gca,'YLim',ylims);
%             set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims)
            ylims = get(gca, 'YLim');
            plot([shift_to shift_to], ylims, 'k', 'LineWidth', 1)
            if vv ==1
                ylabel(sprintf('S%02d',sublist(ss)));
            end

            if ss==1
                if contains(md_names{vv}, ' ')
                    % break it into two strings 
                    spaceind = find(md_names{vv}==' ');
                    title(sprintf('%s\n%s', md_names{vv}(1:spaceind-1), md_names{vv}(spaceind+1:end)));
                else               
                    title(sprintf('%s', md_names{vv}));
                end
            end
            set(gcf,'Color','w')
        end
    end
    % now plot the average across subjects.
    for vv = 1:numel(md_inds)
            
        subplot(nSubj+1,numel(md_inds),(nSubj)*numel(md_inds)+vv);hold all;
      
        for cc = 1:nConds
            recs = squeeze(avg_recs(:,md_inds(vv),cc,:));
            if nSubj>1
                meanvals = nanmean(recs,1);
            else
                meanvals = recs;
            end
            h=errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);
            if vv==numel(vismotor_inds)
                lh = [lh,h];
            end
        end
        
        set(gca, 'FontSize', 12)
         if vv==1
            set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims)
        else
            set(gca,'XTick',[]);
            set(gca,'YLim',ylims);
        end
        ylims = get(gca, 'YLim');
        plot([shift_to shift_to], ylims, 'k', 'LineWidth', 1)
        if vv ==1
            ylabel('Average');           
            ylabel('Avg');
            xlabel('Orientation Channel')    
        end
        
        set(gcf,'Color','w')
    end
    legend(lh,condLabStrs,'Location','SouthEast');
    suptitle('Train Loc, Test WM')
   set(gcf,'Position',[200,200,1400,800]);
end

%% subject average

if plotMDReconsAvg
    lh = [];
    npx = ceil(sqrt(numel(md_inds)));
    npy = ceil(numel(md_inds)/npx);
    figure;hold all;
   
    % now plot the average across subjects.
    for vv = 1:numel(md_inds)
            
        subplot(npx,npy,vv);hold all;
      
        for cc = 1:nConds
            recs = squeeze(avg_recs(:,md_inds(vv),cc,:));
            if nSubj>1
                meanvals = nanmean(recs,1);
            else
                meanvals = recs;
            end
            h= errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);
            if vv==numel(vismotor_inds)
                lh = [lh,h];
            end
        end
        
        set(gca, 'FontSize', fs)
        if vv==numel(md_inds)
            set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims,'YTick',[ylims(1),0,ylims(2)])
            ylabel('a.u.');
            xlabel('Spatial Position (deg)')    
        else
            set(gca,'XLim', [0, 360],'XTick',[]);
            set(gca,'YLim',ylims,'YTick',[ylims(1),0,ylims(2)]);
        end
        title(md_names{vv})
        plot([shift_to shift_to], ylims, 'k', 'LineWidth', 1)
        
        set(gcf,'Color','w')
    end
    legend(lh,condLabStrs,'Location','SouthEast');
    suptitle('Train Loc, Test WM')
    set(gcf,'Position',[200,200,1400,800]);
end

%% 2-way RM anova on fidelity values
% using shuffling to compute significance of each effect
numcores = 8;
nPermIter=1000;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 767657;
[p_vals, ranova_table, iter] = get_f_dist(fidelity(:,vismotor_inds,:), nPermIter, rndseed, 0);

% print results of the shuffling test
f_vals = ranova_table{[3,5,7],[4]}';
df = ranova_table{[3,5,7],[2]}';
array2table([f_vals; df; p_vals],'RowNames',{'f','df','pval'},'VariableNames',{'ROI','Condition','interaction'})

%% now doing pairwise condition comparisons - paired t-test.
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 565646;
rng(rndseed,'twister')

real_sr_stat = nan(nROIs,1);
rand_sr_stat = nan(nROIs, nPermIter);

vals = fidelity;

% p_diff_sr=nan(nROIs,1);
for vv=1:nROIs
    realvals = squeeze(vals(:,vv,:));
    
    % what is the sign-rank statistic for the real data?
    real_sr_stat(vv) = signrank_MMH(realvals(:,1),realvals(:,2));
%     [p,h,stats]=signrank(realvals(:,1),realvals(:,2));
%     p_diff_sr(vv) = p;
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

% compute a two-tailed p-value comparing the real stat to the random
% distribution. Note that the <= and >= are inclusive, because any
% iterations where real==null should count toward the null hypothesis. 
p_diff_sr = 2*min([mean(repmat(real_sr_stat,1,nPermIter)>=rand_sr_stat,2), ...
    mean(repmat(real_sr_stat,1,nPermIter)<=rand_sr_stat,2)],[],2);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

% print out which areas show a significant condition effect across all subj
array2table([diff_is_sig(vismotor_inds), p_diff(vismotor_inds)],'RowNames',vismotor_names,'VariableNames',{'cond_diff','p'})

% print out how many subjects individually showed condition difference in
% the direction of the main effect (random>pred)
array2table(squeeze(sum(vals(:,vismotor_inds,2)>vals(:,vismotor_inds,1),1))','RowNames',vismotor_names,'VariableNames',{'rand_gr_pred'})


 %% plot with single subjects
bw=0.50;
fs=14;

if plotVisMotorFids
   
    vals = fidelity;
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
        b(cc).FaceColor = col(cc,:);
        b(cc).EdgeColor = col(cc,:);
        errorbar(barPos(:,cc),meanVals(:,cc),seVals(:,cc),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    end

    set(gca,'XTick', 1:numel(vismotor_inds))
    set(gca,'XTickLabel', vismotor_names,'XTickLabelRotation',90);
    ylabel('Fidelity')
    set(gca,'YLim',fidlims)
    set(gca,'XLim',[0,numel(vismotor_inds)+1])
   
    set(gca,'FontSize',fs);
    set(gcf,'Position',[800,800,1200,500]);
    % get locations of bars w offsets
    c=get(gcf,'Children');b=get(c(end),'Children');
   
    verspacerbig = range(fidlims)/50;
    horspacer = abs(diff(bar_offset))/2;
%     
    for vv=1:numel(vismotor_inds)
        % add individual subjects
        for ss=1:nSubj
            subvals = squeeze(vals(ss,vismotor_inds(vv),:));
            h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(5,:),'LineWidth',1.5);
            uistack(h,'bottom');
        end
       
        % add significance of condition differences
        for aa=1:numel(alpha_vals)
            if p_diff(vismotor_inds(vv))<alpha_vals(aa)
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
%     uistack(b(end),'top');
%     uistack(b(end-1),'top')
    set(gcf,'color','white')
    set(gcf, 'WindowStyle','normal','WindowState','normal')
    saveas(gcf,fullfile(figpath,'IEM_TrnSWMLoc_allareas.pdf'),'pdf');
end

%% plot fidelity - motor
bw=0.50;
fs=14;

if plotMDFids
   
    vals = fidelity;
    meanvals = squeeze(mean(vals,1));
    semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
    
    meanVals=meanvals(md_inds,:);
    seVals=semvals(md_inds,:);
    
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

    set(gca,'XTick', 1:numel(md_inds))
    set(gca,'XTickLabel', md_names,'XTickLabelRotation',90);
    ylabel('Fidelity')
    set(gca,'YLim',fidlims)
    set(gca,'XLim',[0,numel(md_inds)+1])
   
    set(gca,'FontSize',fs);
    set(gcf,'Position',[800,800,1200,500]);
    % get locations of bars w offsets
    c=get(gcf,'Children');b=get(c(end),'Children');
   
    verspacerbig = range(fidlims)/50;
    horspacer = abs(diff(bar_offset))/2;
%     
    for vv=1:numel(md_inds)
        % add individual subjects
        for ss=1:nSubj
            subvals = squeeze(vals(ss,md_inds(vv),:));
            h=plot(vv+bar_offset,subvals,'.-','Color',sub_colors(5,:),'LineWidth',1.5);
            uistack(h,'bottom');
        end
       
        % add significance of condition differences
        for aa=1:numel(alpha_vals)
            if p_diff(md_inds(vv))<alpha_vals(aa)
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
%     uistack(b(end),'top');
%     uistack(b(end-1),'top')
    set(gcf,'color','white')
    set(gcf, 'WindowStyle','normal','WindowState','normal')
%     saveas(gcf,fullfile(figpath,'IEM_TrainTestWithinConds_allareas.pdf'),'pdf');
end

%%
if plotVisualFidSideMatch
    
    for cc=1:nConds
        vals = squeeze(fidelity_side_match(:,vismotor_inds,cc,:));
        if nSubj>1
            meanvals = squeeze(nanmean(vals,1));
            semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
        else
            meanvals = vals;
            semvals =[];
        end

        plot_barsAndStars(meanvals,semvals,[],[],[],fidlims,vismotor_names,sameSideLabs,'Fidelity',condLabStrs{cc},col_resp)
    end
end


%%
if plotMDFidSideMatch
    
    for cc=1:nConds
        vals = squeeze(fidelity_side_match(:,md_inds,cc,:));
        if nSubj>1
            meanvals = squeeze(nanmean(vals,1));
            semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
        else
            meanvals = vals;
            semvals =[];
        end

        plot_barsAndStars(meanvals,semvals,[],[],[],fidlims,md_names,sameSideLabs,'Fidelity',condLabStrs{cc},col_resp)
    end
end



%%
if plotVisBiasSepDirs

    npx = ceil(sqrt(numel(vismotor_inds)));
    npy = ceil(numel(vismotor_inds)/npx);
    ylims_small=[-0.5, 0.5];
    dir_colors = plasma(nDirs+2);

    for cc=1:nConds
        for bb=1:nBins
            lh=[];
            figure;hold all;    
            for vv=1:numel(vismotor_inds)

                subplot(npx,npy,vv);hold all;

                for ww=1:nDirs
                    recs = squeeze(recs_binned_preview_bound(:,vismotor_inds(vv),cc,bb,ww,:));
                    if nSubj>1
                        meanvals = nanmean(recs,1);
                    else
                        meanvals = recs;
                    end
                    [~,peakind] = max(meanvals);
                    h= errorbar(xx',meanvals,[], 'Color',dir_colors(ww,:),'LineWidth',1);
                    plot([xx(peakind),xx(peakind)], ylims_small,'Color',dir_colors(ww,:));
                    if vv==numel(vismotor_inds)
                        lh = [lh,h];
                    end
                end

                set(gca, 'FontSize', 12)
                if vv==1
                    set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims_small)
                else
                    set(gca,'XTick',[]);
                    set(gca,'YLim',ylims_small);
                end

                plot([shift_to shift_to], ylims_small, 'k', 'LineWidth', 1)
                if vv ==1
            %         ylabel('Average');           
                    ylabel('Avg');
                    xlabel('Orientation Channel')    
                end

                if contains(vismotor_names{vv}, ' ')
                    % break it into two strings 
                    spaceind = find(vismotor_names{vv}==' ');
                    title(sprintf('%s\n%s', vismotor_names{vv}(1:spaceind-1), vismotor_names{vv}(spaceind+1:end)));
                else               
                    title(sprintf('%s', vismotor_names{vv}));
                end
                set(gcf,'Color','w')
            end

            legend(lh,dirLabs,'Location','SouthEast');
            suptitle(sprintf('%s\nTarget %s',condLabStrs{cc},binLabs{bb}))
            set(gcf,'Position',[200,200,1400,800]);

        end
    end
end
%%
if plotVisBiasAvgDirs
    
npx = ceil(sqrt(numel(vismotor_inds)));
npy = ceil(numel(vismotor_inds)/npx);
ylims_small=[-0.5, 0.5];
bin_colors = plasma(nBins);

for bb=1:nBins
    lh=[];
    figure;hold all;    
    for vv=1:numel(vismotor_inds)

        subplot(npx,npy,vv);hold all;

        for cc=1:nConds
            recs1 = squeeze(recs_binned_preview_bound(:,vismotor_inds(vv),cc,bb,1,:));
            recs2 = squeeze(recs_binned_preview_bound(:,vismotor_inds(vv),cc,bb,2,:));
            recs = mean(cat(3,flipud(recs1), recs2),3);
            if nSubj>1
                meanvals = nanmean(recs,1);
            else
                meanvals = recs;
            end
            [~,peakind] = max(meanvals);
            h= errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);
            plot([xx(peakind),xx(peakind)], ylims_small,'Color',col(cc,:));
            if vv==numel(vismotor_inds)
                lh = [lh,h];
            end
        end

        set(gca, 'FontSize', 12)
        if vv==1
            set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims_small)
        else
            set(gca,'XTick',[]);
            set(gca,'YLim',ylims_small);
        end

        plot([shift_to shift_to], ylims_small, 'k', 'LineWidth', 1)
        if vv ==1
    %         ylabel('Average');           
            ylabel('Avg');
            xlabel('Orientation Channel')    
        end
        if contains(vismotor_names{vv}, ' ')
            % break it into two strings 
            spaceind = find(vismotor_names{vv}==' ');
            title(sprintf('%s\n%s', vismotor_names{vv}(1:spaceind-1), vismotor_names{vv}(spaceind+1:end)));
        else               
            title(sprintf('%s', vismotor_names{vv}));
        end
        set(gcf,'Color','w')
    end

    legend(lh,condLabStrs,'Location','SouthEast');
    suptitle(sprintf('%s\nRight is toward category center',binLabs{bb}))
    set(gcf,'Position',[200,200,1400,800]);

end

end

%%
if plotVisualBiasBars

    for bb=1:nBins
        vals = squeeze(bias(:,vismotor_inds,:,bb));       
        if nSubj>1
            meanvals = squeeze(nanmean(vals,1));
            semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
        else
            meanvals = squeeze(vals);
            semvals =[];
        end
        plot_barsAndStars(meanvals,semvals,[],[],[],biaslims,vismotor_names,condLabStrs,'Avg bias (flipped so up is away)',sprintf('Bias away from bound\n%s',binLabs{bb}),col(1:2,:))
    end
end

%%
if plotMotorBiasBars
    
    for bb=1:nBins
        vals = squeeze(bias(:,md_inds,:,bb));       
        if nSubj>1
            meanvals = squeeze(nanmean(vals,1));
            semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
        else
            meanvals = squeeze(vals);
            semvals =[];
        end
        plot_barsAndStars(meanvals,semvals,[],[],[],biaslims,md_names,condLabStrs,'Avg bias (flipped so up is away)',sprintf('Bias away from bound\n%s',binLabs{bb}),col(1:2,:))

    end
end