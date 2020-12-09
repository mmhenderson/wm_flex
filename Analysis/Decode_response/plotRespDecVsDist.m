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
    'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

nROIs = length(ROI_names);

plot_order1 = [1:5,10,11,6:9];  % visual ROIs
plot_order2 = [12:14];   % motor ROIs
plot_order3 = [15:20];

visual_names = ROI_names(plot_order1);
motor_names = ROI_names(plot_order2);
md_names = ROI_names(plot_order3);

plot_order_all = [plot_order1, plot_order2, plot_order3];

vis_inds = find(ismember(plot_order_all,plot_order1));
motor_inds = find(ismember(plot_order_all,plot_order2));
md_inds = find(ismember(plot_order_all,plot_order3));

nVox2Use = 10000;
nPermIter=1000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 0.9];
conflims=[-0.5, 2.5];
dprimelims = [-0.2, 1.4];



condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

chance_val=0.5;

plotMotorConf = 1;
plotMotorAcc=1;
plotMotorAccOverTime=1;
plotMotorConfOverTime=1;
plotMotorScatterConf=1;

% when looking at distance from boundary effects, how many distance bins to use?
nDistBins=3;
dist_bin_edges = linspace(0,90,nDistBins+1);
distBinLabs=[];
dist_bin_centers = dist_bin_edges(1:end-1)+diff(dist_bin_edges(1:2))/2;
for bb=1:nDistBins
    distBinLabs{bb} = sprintf('%d-%d deg from bound', dist_bin_edges(bb),dist_bin_edges(bb+1));
end
if nDistBins==1
    distBinLabs={''};
end
nDirs =  2; % CW/CCW
dirLabs = {'Targ CW (-) of bound','Targ CCW (+) of bound'};

nSpatBins=8;
spat_bin_centers=[0:45:359];
spat_bin_size=diff(spat_bin_centers(1:2));

col_conds = plasma(5);
col_conds = col_conds(2:2:end-1,:);
if nDistBins>1
    scale_vec = linspace(-.2,.2,nDistBins);
    scale_mat = permute(repmat(scale_vec,2,1,3),[1,3,2]);
    col_conds = max(min(scale_mat + repmat(col_conds,1,1,numel(scale_vec)),1),0);
end

diff_col=[0.5, 0.5, 0.5];
ms=10;  % marker size for significance dots
nTrialsTotal = 400;
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
% events to plot as vertical lines
evts2plot = [3.5, 4.5, 16.5, 18.5];
chance_val=0.5;
lw=1;

fs=12;
acclims=[0.4, 1];
%% load results

confbins_allsubs = nan(nSubj,nROIs,nConds,nDistBins);
accbins_allsubs = nan(nSubj,nROIs,nConds,nDistBins);
accbinstime_allsubs = nan(nSubj,nROIs,nConds,nDistBins,nTRs_out);
confbinstime_allsubs = nan(nSubj,nROIs,nConds,nDistBins,nTRs_out);
confall = nan(nSubj, nROIs, nTrialsTotal);
condall=nan(nSubj, nTrialsTotal);
distbinall = nan(nSubj, nTrialsTotal);
distall = nan(nSubj, nTrialsTotal);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
%     fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
%     load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    allconf1=allconf;
    
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    % get some labels about all the trials 
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    % how far each trial is from boundary, binned
    dist_to_bound = mainSig(1).dist_to_real_bound;
    dist_to_bound(condlabs==2) = mainSig(1).dist_to_rand_bound(condlabs==2);
    which_bin = zeros(size(dist_to_bound));
    for bb=1:nDistBins
        inds2use = (dist_to_bound<dist_bin_edges(bb+1) &...
                   dist_to_bound>dist_bin_edges(bb));
        which_bin(inds2use) = bb;
    end
    condlabs = mainSig(1).condLabs;
    
    confall(ss,:,:) = allconf1;
    condall(ss,:) = condlabs;
    distbinall(ss,:) = dist_to_bound;
    distall(ss,:) = dist_to_bound;
    % real trial labels, bins of spatial position
%     tstPosLabs = mainSig(vv).targPos;
%     tstBinLabs = zeros(size(tstPosLabs));        
%     for bb=1:nSpatBins
%         inds_this_bin = abs(tstPosLabs-(spat_bin_centers(bb)-0.0001))<spat_bin_size/2 | abs((tstPosLabs-360)-(spat_bin_centers(bb)-0.0001))<spat_bin_size/2;
%         tstBinLabs(inds_this_bin) = bb;
%     end
%     assert(~any(tstBinLabs==0))
%     neach_tst = sum(repmat(tstBinLabs,1,nSpatBins)==repmat((1:nSpatBins),size(tstBinLabs,1),1));
    
    % calculate performance within particular bins of target-bbound
    % distance.
    for bb=1:nDistBins
        for cc=1:nConds
            cvals = allconf1(:,condlabs==cc & which_bin==bb);
            confbins_allsubs(ss,:,cc,bb) = mean(cvals,2);   % average over trials
            accbins_allsubs(ss,:,cc,bb) = mean(cvals>0,2);
            for tr=1:nTRs_out
                cvals = allconf(:,condlabs==cc & which_bin==bb,tr);
                confbinstime_allsubs(ss,:,cc,bb,tr) = mean(cvals,2);
                accbinstime_allsubs(ss,:,cc,bb,tr) = mean(cvals>0,2);
            end
        end
    end
   
end

assert(~any(isnan(confbins_allsubs(:))))
assert(~any(isnan(accbins_allsubs(:))))

%% make a bar plot of acc - visual areas
if plotMotorConf
    
   for cc=1:nConds
       vals = squeeze(confbins_allsubs(:,motor_inds,cc,:));       
       meanvals = squeeze(mean(vals,1));
       semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
        
       plot_barsAndStars(meanvals,semvals,[],...
        [],0,conflims,motor_names,distBinLabs,...
        'Classifier Confidence',sprintf('Correct Resp\n%s',condLabStrs{cc}),squeeze(col_conds(cc,:,:))')
       set(gcf,'Position',[800,800,1200,420])
   end

end

%%
if plotMotorScatterConf
     
    sub_colors = plasma(nSubj+1);
    conflims=[-12,12];

    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(numel(motor_inds)));
        npy = ceil(numel(motor_inds)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(motor_inds)
            subplot(npx,npy,vv);hold all;

%             for ss=1
            for ss=1:nSubj

                x=squeeze(distall(ss,condall(ss,:)==cc));
                y=squeeze(confall(ss,motor_inds(vv),condall(ss,:)==cc));
                
                plot(x,y,'.','Color',sub_colors(ss,:));
            end
           
            set(gca, 'FontSize', 12)
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
            set(gca,'XLim',[0,90],'YLim',conflims);
            plot([0,90],[0,0],'-','Color','k');
            plot([0,0],conflims,'-','Color','k');
            axis square
            xlabel('Distance from bound (deg)')
            ylabel('Response decoder confidence')
            
            title(sprintf('%s', motor_names{vv}));
            
        end

        suptitle(sprintf('Resp decoding in each area vs. target distance from bound\n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end
%%
if plotMotorAcc
    
   for cc=1:nConds
       vals = squeeze(accbins_allsubs(:,motor_inds,cc,:));       
       meanvals = squeeze(mean(vals,1));
       semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
        
       plot_barsAndStars(meanvals,semvals,[],...
        [],chance_val,acclims,motor_names,distBinLabs,...
        'Accuracy',sprintf('Correct Resp\n%s',condLabStrs{cc}),squeeze(col_conds(cc,:,:))')
       set(gcf,'Position',[800,800,1200,420])
   end

end

%%
if plotMotorAccOverTime
    for cc = 1:nConds
        nplots_motor = ceil(sqrt(numel(motor_inds)));
        figure();hold all;

        for vi = 1:numel(motor_inds)

            subplot(nplots_motor,ceil(numel(motor_inds)/nplots_motor),vi);hold all;

            lh=[];
            for bb=1:nDistBins
%             meanvals = meanvals_allcond(cc,:);
%             semvals = sem_within_allcond(cc,:);
                vals = squeeze(accbinstime_allsubs(:,motor_inds(vi),cc,bb,:));
                if nSubj==1
                    meanvals =squeeze(vals)';
                    semvals = [];
                else
                    meanvals = squeeze(mean(vals,1))';
                    semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
                end

                lh=[lh,plot(tax,meanvals,'-','Color',squeeze(col_conds(cc,:,bb))','LineWidth',lw)];

                bandedError_MMH(tax, meanvals',semvals',squeeze(col_conds(cc,:,bb)), 0.5);
    %             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
     
            end
        
            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(acclims);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Accuracy');
            end

            if vi==numel(motor_inds)
                legend(lh,distBinLabs,'FontSize', fs);
            end

            title(sprintf('%s', motor_names{vi}));
        end
      
    
        suptitle(condLabStrs{cc});
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
    end
end

%%
if plotMotorConfOverTime
    for cc = 1:nConds
        nplots_motor = ceil(sqrt(numel(motor_inds)));
        figure();hold all;

        for vi = 1:numel(motor_inds)

            subplot(nplots_motor,ceil(numel(motor_inds)/nplots_motor),vi);hold all;

            lh=[];
            for bb=1:nDistBins
%             meanvals = meanvals_allcond(cc,:);
%             semvals = sem_within_allcond(cc,:);
                vals = squeeze(confbinstime_allsubs(:,motor_inds(vi),cc,bb,:));
                if nSubj==1
                    meanvals =squeeze(vals)';
                    semvals = [];
                else
                    meanvals = squeeze(mean(vals,1))';
                    semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
                end

                lh=[lh,plot(tax,meanvals,'-','Color',squeeze(col_conds(cc,:,bb))','LineWidth',lw)];

                bandedError_MMH(tax, meanvals',semvals',squeeze(col_conds(cc,:,bb)), 0.5);
    %             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
     
            end
        
            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(conflims);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],conflims,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Classifier confidence (a.u.)');
            end

            if vi==numel(vis_inds)
                legend(lh,distBinLabs,'FontSize', fs);
            end

            title(sprintf('%s', motor_names{vi}));
        end
      
    
        suptitle(condLabStrs{cc});
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
    end
end
