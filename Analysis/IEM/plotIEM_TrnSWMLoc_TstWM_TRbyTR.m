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

plotVisMotorFids = 1;
plotVisMotorFidSideMatch = 0;
plotMDFids = 0;
plotMDFidSideMatch=0;

plotVisualBias=0;
plotMotorMDBias=0;
plotVisBiasSepDirs=0;
plotVisRecCenterBound=0;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

nVOIs = length(plot_order_all);

trDur = 0.80;

% col = plasma(5);
% col = col(2:2:end-1,:);
col = [125, 93, 175; 15, 127, 98]./255;

ylims_fid = [-0.1, 0.25];
ylims_bias=[-90,90];
lw = 1;
fs=20;

alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,14,20];
alpha=alpha_vals(1);


col_resp = plasma(3);
col_resp = col_resp(1:2,:);
% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.19, 0.21, 0.23];
diff_col=[0.5, 0.5, 0.5];



alpha_vals=[0.05,0.01,0.001];
alpha_ms = [8,16,24];
% when looking at category repulsion effects, how many distance bins to use?
nBins=1;
dist_bin_edges = linspace(0,90,nBins+1);
for bb=1:nBins
    binLabs{bb} = sprintf('%d-%d deg from bound', dist_bin_edges(bb),dist_bin_edges(bb+1));
end
nDirs =  2; % CW/CCW
dirLabs = {'Bound CCW (+) of target','Bound CW (-) of target'};

sameSideLabs = {'Finger on same side as target','Finger on opp side of target'};

%%
for ss=1:nSubj
    
    substr = sprintf('S%02d',sublist(ss));    
    fn = fullfile(exp_path,'Analysis','IEM','IEM_results',sprintf('TrnSWMLoc_TestWM_TRbyTR_%s.mat',substr));
    fprintf('loading from %s\n',fn);
    load(fn);
    assert(numel(allchanresp)==numel(ROI_names));
    if ss==1        
        nTRs_out = size(allchanresp(1).chan_resp,3);
        tax = trDur*(0:nTRs_out-1);
        %preallocate some arrays       
        fidelity = nan(nSubj, nVOIs, nConds, nTRs_out);
        fidelity_side_match = nan(nSubj, nVOIs, nConds, nTRs_out, length(sameSideLabs));
        bias = nan(nSubj, nVOIs, nConds, nBins, nTRs_out);
        recs_binned_preview_bound = nan(nSubj, nVOIs, nConds, nBins, nDirs, length(xx), nTRs_out);
        recs_centered_preview_bound = nan(nSubj, nVOIs, nConds, length(xx), nTRs_out);
        
    end
    
%     respLabs=allchanresp(1).respLabs;
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

    % check that the above code is working by making sure it gets the
    % right answer for the predictable trials
%     assert(all(resp_assoc_preview(condlabs==1)==respLabs(condlabs==1)));

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
           
           for tr = 1:nTRs_out
           
               % take out trials from one condition at a time
               theserecs = allchanresp(plot_order_all(vv)).chan_resp_shift(condlabs==cc,:,tr);
               meanrec = squeeze(mean(theserecs,1));

               % get the fidelity
               angs = abs((1:360)-180);
               cos_vals = cosd(angs);
               fidelity(ss,vv,cc,tr) = mean(cos_vals.*meanrec);

               if ss>1 && fidelity(ss-1,vv,cc,tr)==fidelity(ss,vv,cc,tr)
                   error('you have a missing value and something bad is happening')
               end
               
               for bb=1:nBins
                   
                   for ww=1:nDirs
               
                       % ww=2 for CCW, ww=1 for CW
                       inds2use = which_bin(condlabs==cc)==bb & which_dir(condlabs==cc)==ww;
                       recs_this_bin = theserecs(inds2use,:);
                       recs_binned_preview_bound(ss,vv,cc,bb,ww,:,tr) = mean(recs_this_bin,1);

                   end
                   % ww=1 for CW, ww=2 for CCW
                   indsCW = which_bin(condlabs==cc)==bb & which_dir(condlabs==cc)==1;
                   indsCCW = which_bin(condlabs==cc)==bb & which_dir(condlabs==cc)==2;
                   recsCW = squeeze(mean(theserecs(indsCW,:)))';
                   recsCCW = squeeze(mean(theserecs(indsCCW,:)))';
 
                   rec = mean(cat(2,flipud(recsCW), recsCCW),2);
                   [~,peakind] = max(rec);
                   bias(ss,vv,cc,bb,tr) = xx(peakind)-180;
                   
               end
               
               % now compute recons centered on the preview disk orientation.
               % where do you want them to be centered? this works well b/c center of
               % xx axis.
               shift_to = 180;

               theserecs_unshifted = allchanresp(plot_order_all(vv)).chan_resp(condlabs==cc,:,tr);
               shift_to_bound = zeros(size(theserecs_unshifted));        
               % the final location to shift them to is the preview
               % boundary location (0-360), with 180 degrees added if the
               % response is 2. this makes it so that the half of the disk
               % associated with the response is always CCW (+) of the
               % boundary position.
               bound_pos_to_shift = preview_bound_pos(condlabs==cc);
               bound_pos_to_shift(resp_assoc_preview(condlabs==cc)==2) = bound_pos_to_shift(resp_assoc_preview(condlabs==cc)==2) + 180;
               for tt = 1:size(theserecs_unshifted,1)
                
                    rec = squeeze(theserecs_unshifted(tt,:));
                    shift_to_bound(tt,:) =  wshift('1D', rec,ceil(bound_pos_to_shift(tt)-shift_to));
                
               end
               recs_centered_preview_bound(ss,vv,cc,:,tr)  =mean(shift_to_bound,1);
               
               
               % compute fidelity separately for trials where target was on
               % same side as response (e.g. left finger/left side of screen)
               for mm=1:2
                   theserecs = squeeze(allchanresp(plot_order_all(vv)).chan_resp_shift(condlabs==cc & same_side==mm,:,tr));           
                   fidelity_side_match(ss,vv,cc,tr,mm) = mean(cos_vals.*mean(theserecs,1));
               end
               
           end
           
           
       end
    end
  
end

% compute some basic stats here
vals = fidelity;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
%% now doing pairwise condition comparisons - paired t-test.
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 143545;
rng(rndseed,'twister')
nPermIter=1000;
real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);

% p_diff_sr=nan(nROIs,nTRs_out);
for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,tr));
%         [p,h,stats]=signrank(realvals(:,1),realvals(:,2));
%         p_diff_sr(vv,tr) = p;
        % what is the sign-rank statistic for the real data?
        real_sr_stat(vv,tr) = signrank_MMH(realvals(:,1),realvals(:,2));

        % determine before the parfor loop which conditions get randomly
        % swapped on each iteration (otherwise not deterministic)
        inds2swap = double(randn(nSubj,nPermIter)>0);
        inds2swap(inds2swap==0) = -1;

        parfor ii=1:nPermIter          

            % randomly permute the condition labels within subject
            randvals=realvals;
            randvals(inds2swap(:,ii)==-1,:) = randvals(inds2swap(:,ii)==-1,[2,1]);    
            % what is the sign-rank statistic for this randomly permuted data?
            rand_sr_stat(vv,tr,ii) = signrank_MMH(randvals(:,1),randvals(:,2));

        end
    end
end

% % compute a two-tailed p-value comparing the real stat to the random
% % distribution. Note that the <= and >= are inclusive, because any
% % iterations where real==null should count toward the null hypothesis. 
p_diff_sr = 2*min(cat(3,mean(repmat(real_sr_stat,1,1,nPermIter)>=rand_sr_stat,3), ...
    mean(repmat(real_sr_stat,1,1,nPermIter)<=rand_sr_stat,3)),[],3);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

array2table(squeeze(p_diff_sr(vismotor_inds,1:15)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Diff_TR',1:15))
array2table(squeeze(p_diff_sr(vismotor_inds,16:30)),...
    'RowNames',vismotor_names,'VariableNames',strseq('Diff_TR',16:30))

%% Now we make a figure for all subjects
if plotVisMotorFids
    
    nplots_vis = ceil(sqrt(numel(vismotor_inds)));
    figure();hold all;
    
    for vi = 1:numel(vismotor_inds)
       
        subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;
        
        lh=[];
        for cc = 1:nConds

            valsplot = fidelity(:,vismotor_inds(vi),cc,:);
            if nSubj==1
                meanvals =squeeze(valsplot)';
                semvals = [];
            else
                meanvals = squeeze(mean(valsplot,1))';
                semvals = squeeze(std(valsplot,[],1))'./sqrt(nSubj);
            end
            lh=[lh,plot(tax,meanvals,'-','Color',col(cc,:),'LineWidth',lw)];
            bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);

        end
        ll=condLabStrs;
        for aa=1:numel(alpha_vals)
            inds2plot=p_diff(vismotor_inds(vi),:)<alpha_vals(aa);
            plot(tax(inds2plot), repmat(sig_heights(nConds+1),sum(inds2plot),1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))
            lh=[lh, plot(tax(1)-5, repmat(sig_heights(nConds+1),1,1),'.','Color',diff_col,'MarkerSize',alpha_ms(aa))];
            ll{numel(ll)+1} = sprintf('p<%.03f',alpha_vals(aa));
        end
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(ylims_fid);
        h=plot(get(gca,'XLim'),[0,0],'--','Color',[0.7, 0.7, 0.7]);
        uistack(h,'bottom')
        for ee = 1:length(evts2plot)
            h=plot([evts2plot(ee),evts2plot(ee)],ylims_fid,'-','Color',[0.8, 0.8, 0.8]);
            uistack(h,'bottom')
        end
        h=fill([evts2plot(1),evts2plot(2), evts2plot(2),evts2plot(1)], repelem(ylims_fid,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        h=fill([evts2plot(3),evts2plot(4), evts2plot(4),evts2plot(3)], repelem(ylims_fid,2),[0.8, 0.8, 0.8],'EdgeColor','none');
        uistack(h,'bottom')
        
        if vi==1
            xlabel('Time(s)');
            ylabel('Fidelity');
        end

        if vi==numel(vismotor_inds)
            legend(lh,ll,'FontSize', fs);
        end
       
        title(sprintf('%s', vismotor_names{vi}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200]);
end

%%
if plotVisMotorFidSideMatch
    
    for cc=1:nConds
        figure;hold all;
        nplots_vis = ceil(sqrt(numel(vismotor_inds)));

        for vi = 1:numel(vismotor_inds)
            subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;


            for mm=1:2
                vals = fidelity_side_match(:,vismotor_inds(vi),cc,:,mm);
                if nSubj==1
                    meanvals =squeeze(vals);
                    semvals = [];
                else
                    meanvals = squeeze(nanmean(vals,1));
                    semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
                end

%                 errorbar(tax, meanvals, semvals,'Color',col_resp(mm,:),'LineWidth',lw);
                bandedError_MMH(tax, meanvals,semvals, col_resp(mm,:), 0.5);
            end

            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(ylims_fid);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],ylims_fid,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==numel(vismotor_inds)-1
                xlabel('Time(s)');
                ylabel('Fidelity');
            end

            if vi==numel(vismotor_inds)
                legend(sameSideLabs);
            end

            if contains(vismotor_names{vi}, ' ')
                % break it into two strings
                spaceind = find(vismotor_names{vi}==' ');
                title(sprintf('%s\n%s', vismotor_names{vi}(1:spaceind-1), vismotor_names{vi}(spaceind+1:end)));
            else
                title(sprintf('%s', vismotor_names{vi}));
            end


        end
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,1000]);
        suptitle(condLabStrs{cc});
    end
end
%%
if plotVisualBias
    
    for bb=1:nBins
        figure;hold all;
        nplots_vis = ceil(sqrt(numel(vismotor_inds)));

        for vi = 1:numel(vismotor_inds)
            subplot(nplots_vis,ceil(numel(vismotor_inds)/nplots_vis),vi);hold all;


            for cc=1:nConds
                vals = bias(:,vismotor_inds(vi),cc,bb,:);
                if nSubj==1
                    meanvals =squeeze(vals);
                    semvals = [];
                else
                    meanvals = squeeze(nanmean(vals,1));
                    semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
                end

%                 errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
                bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);
            end

            set(gca, 'FontSize', 20, 'XLim',[0 max(tax)])
            ylim(ylims_bias);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],ylims_bias,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Bias (+ is away from bound)');
            end

            if vi==numel(vismotor_inds)
                legend(condLabStrs);
            end

            if contains(vismotor_names{vi}, ' ')
                % break it into two strings
                spaceind = find(vismotor_names{vi}==' ');
                title(sprintf('%s\n%s', vismotor_names{vi}(1:spaceind-1), vismotor_names{vi}(spaceind+1:end)));
            else
                title(sprintf('%s', vismotor_names{vi}));
            end


        end
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
        suptitle(binLabs{bb});
    end
end
%%
if plotMDFids
    
    figure;hold all;
    nplots_motor = ceil(sqrt(numel(md_inds)));
        
    for vi = 1:numel(md_inds)
        subplot(nplots_motor,ceil(numel(md_inds)/nplots_motor),vi);hold all;
        
        lh=[];
        for cc=1:nConds
            vals = fidelity(:,md_inds(vi),cc,:);
            if nSubj==1
                meanvals =squeeze(vals)';
                semvals = [];
            else
                meanvals = squeeze(nanmean(vals,1))';
                semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)))';
            end
            lh=[lh,plot(tax,meanvals,'Color',col(cc,:),'LineWidth',lw)];
%             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
            bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);
        end
        
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(ylims_fid);
        plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],ylims_fid,'-','Color',[0.8, 0.8, 0.8]);
        end
        
        if vi==numel(vismotor_inds)-1
            xlabel('Time(s)');
            ylabel('Fidelity');
        end

        if vi==numel(md_inds)
            legend(lh,condLabStrs);
        end
        
        if contains(md_names{vi}, ' ')
            % break it into two strings
            spaceind = find(md_names{vi}==' ');
            title(sprintf('%s\n%s', md_names{vi}(1:spaceind-1), md_names{vi}(spaceind+1:end)));
        else
            title(sprintf('%s', md_names{vi}));
        end
        
        
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,1000]);
end

%%
if plotMDFidSideMatch
    
    for cc=1:nConds
        figure;hold all;
        nplots_motor = ceil(sqrt(numel(md_inds)));

        for vi = 1:numel(md_inds)
            subplot(nplots_motor,ceil(numel(md_inds)/nplots_motor),vi);hold all;


            for mm=1:2
                vals = fidelity_side_match(:,md_inds(vi),cc,:,mm);
                if nSubj==1
                    meanvals =squeeze(vals);
                    semvals = [];
                else
                    meanvals = squeeze(nanmean(vals,1));
                    semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
                end

%                 errorbar(tax, meanvals, semvals,'Color',col_resp(mm,:),'LineWidth',lw);
                bandedError_MMH(tax, meanvals,semvals, col_resp(mm,:), 0.5);
            end

            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(ylims_fid);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],ylims_fid,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==numel(vismotor_inds)-1
                xlabel('Time(s)');
                ylabel('Fidelity');
            end

            if vi==numel(md_inds)
                legend(sameSideLabs);
            end

            if contains(md_names{vi}, ' ')
                % break it into two strings
                spaceind = find(md_names{vi}==' ');
                title(sprintf('%s\n%s', md_names{vi}(1:spaceind-1), md_names{vi}(spaceind+1:end)));
            else
                title(sprintf('%s', md_names{vi}));
            end


        end
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,1000]);
        suptitle(condLabStrs{cc});
    end
end

%%
if plotMotorMDBias
    
    for bb=1:nBins
        figure;hold all;
        nplots_motor = ceil(sqrt(numel(md_inds)));

        for vi = 1:numel(md_inds)
            subplot(nplots_motor,ceil(numel(md_inds)/nplots_motor),vi);hold all;


            for cc=1:nConds
                vals = bias(:,md_inds(vi),cc,bb,:);
                if nSubj==1
                    meanvals =squeeze(vals);
                    semvals = [];
                else
                    meanvals = squeeze(nanmean(vals,1));
                    semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
                end

%                 errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
                bandedError_MMH(tax, meanvals,semvals, col(cc,:), 0.5);
            end

            set(gca, 'FontSize', 12, 'XLim',[0 max(tax)])
            ylim(ylims_bias);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],ylims_bias,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Bias (+ is away from bound)');
            end

            if vi==numel(md_inds)
                legend(condLabStrs);
            end

            if contains(md_names{vi}, ' ')
                % break it into two strings
                spaceind = find(md_names{vi}==' ');
                title(sprintf('%s\n%s', md_names{vi}(1:spaceind-1), md_names{vi}(spaceind+1:end)));
            else
                title(sprintf('%s', md_names{vi}));
            end


        end
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
        suptitle(binLabs{bb});
    end
end

%%
tr2plot = 10;
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
                    recs = squeeze(recs_binned_preview_bound(:,vismotor_inds(vv),cc,bb,ww,:,tr2plot));
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
                if vv==numel(vismotor_inds)-1
                    set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims_small)
                    ylabel('Avg');
                    xlabel('Orientation Channel')    
                else
                    set(gca,'XTick',[]);
                    set(gca,'YLim',ylims_small,'XLim', [0, 360]);
                end

                plot([shift_to shift_to], ylims_small, 'k', 'LineWidth', 1)
%                 plot([270,270], ylims_small, 'Color',[0.8, 0.8, 0.8], 'LineWidth', 1)


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
            suptitle(sprintf('%s',condLabStrs{cc}))
            set(gcf,'Position',[200,200,1400,800]);

        end
    end
end
%%
tr2plot = 10;
if plotVisRecCenterBound
%     for tr2plot=[8:20]
    for tr2plot=[10]
        npx = ceil(sqrt(numel(vismotor_inds)));
        npy = ceil(numel(vismotor_inds)/npx);
        ylims_small=[-0.5, 0.5];
    %     dir_colors = plasma(nDirs+2);



        lh=[];
        figure;hold all;    
        for vv=1:numel(vismotor_inds)

            subplot(npx,npy,vv);hold all;

            for cc=1:nConds
                recs = squeeze(recs_centered_preview_bound(:,vismotor_inds(vv),cc,:,tr2plot));
                if nSubj>1
                    meanvals = nanmean(recs,1);
                else
                    meanvals = recs;
                end
                h= errorbar(xx',meanvals,[], 'Color',col(cc,:),'LineWidth',1);

            end 


            set(gca, 'FontSize', 12)
            if vv==numel(vismotor_inds)-1
                set(gca, 'XLim', [0, 360],'XTick',[0:180:360],'XTickLabel',{-180:180:180}, 'YLim',ylims_small)
                ylabel('Avg');
                xlabel('Orientation Channel')    
            else
                set(gca,'XTick',[]);
                set(gca,'YLim',ylims_small,'XLim', [0, 360]);
            end

            plot([shift_to shift_to], ylims_small, 'k', 'LineWidth', 1)
            plot([270,270], ylims_small, 'Color',[0.8, 0.8, 0.8], 'LineWidth', 1)

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
        suptitle(sprintf('Shifted to preview bound position\n%.1f sec after targ onset',tax(tr2plot)))
        set(gcf,'Position',[200,200,1400,800]);
    end
end
    