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

spat_areas = [1:5,10,11,6:9];  % visual ROIs
% resp_areas = [12:15,22];

nROIs_spat = length(spat_areas);
% nROIs_resp = length(resp_areas);

nVox2Use = 10000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

plotSpatAccVsAllAccBehav=0;
plotSpatAccVsRunAccBehav=0;

plotSpatAccVsAllRTBehav=1;
plotSpatAccVsRunRTBehav=1;
plotSpatAccVsRunRTBehavZscore=1;
plotSpatAccVsRunRTBoot=1;

plotSpatConfVsRT = 0;
plotSpatConfVsRTZscore=1;
plotSpatConfVsRTBins = 1;
plotSpatConfVsCorrect = 0;
%% load results
nTrialsTotal = 400;
nRuns=20;
nTrialsPerRun=20;


spatacc_allsubs = nan(nSubj,nROIs_spat,nConds);
% respacc_allsubs = nan(nSubj,nROIs_resp,nConds);
spatconf = nan(nSubj, nROIs_spat, nTrialsTotal);
% respconf = nan(nSubj, nROIs_resp, nTrialsTotal);
rt_allsubs = nan(nSubj, nTrialsTotal);
behavcorrect_allsubs = nan(nSubj, nTrialsTotal);
condlabs_allsubs = nan(nSubj, nTrialsTotal);

runlabs=repelem(1:nRuns, nTrialsPerRun);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    % loading performance of spatial position classifier
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    fn2load = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
%     fn2load = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_TRbyTR_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
%     load(fn2load);
    
    spatacc_allsubs(ss,:,:) = mean(squeeze(allacc(spat_areas,:,:)),3);
    spatconf(ss,:,:) = allconf(spat_areas,:);
    rt_allsubs(ss,:) = rt;
    behavcorrect_allsubs(ss,:) = correct;
    condlabs_allsubs(ss,:) = condlabs;
    clear allacc
    clear allconf
    
    % loading performance of spatial classifier
%     save_dir = fullfile(exp_path,'Analysis','Decode_spatial','Decoding_results');
%     fn2load = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
%     load(fn2load);
%     
%     respacc_allsubs(ss,:,:) = allacc(resp_areas,:);
%     respconf(ss,:,:) = allconf(resp_areas,:);
%     clear allacc
%     clear allconf
end

assert(~any(isnan(spatacc_allsubs(:))))
% assert(~any(isnan(respacc_allsubs(:))))

%% checking to make sure the confidence values give the correct accuracy calculations
% note they aren't perfectly identical (but should be very close) due to
% the way that I averaged over groups to get accuracy and then averaged
% those accuracy values to get teh final accuracy (what is stored in allacc
% variable) 
acc_check = nan(nSubj, numel(spat_areas), nConds);
for ss=1:nSubj
    for vv=1:numel(spat_areas)
        for cc=1:nConds
            conf = spatconf(ss,vv,condlabs_allsubs(ss,:)==cc);
            acc_check(ss,vv,cc) = mean(conf>0);
        end
    end
end
vals = acc_check;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));


%% bootstrap correlations of RT versus spat decoding performance (using run-averaged values)
nPermIter=1000;
shuff_corrs = nan(nROIs_spat, nConds, nPermIter);
real_corrs = nan(nROIs_spat, nConds);

for vv = 1:nROIs_spat
    for cc=1:nConds
        % for this ROI and cond, want to get a single value for each
        % subject and each run of each task (e.g. 20*6 total)
        allx=[];ally=[];
        slabs=[];
        for run=1:nRuns
            means=zeros(nSubj, 2);
            for ss=1:nSubj
                rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                
                % estimate accuracy on trials in this run
                x = mean(sconf>0);
                y = median(rts);
                allx=[allx;x];
                ally=[ally;y];
                slabs=[slabs;ss];
            end
        end
        
     
       rho = corr(allx,ally);
       real_corrs(vv,cc) = rho;
       % now same thing for shuffled data
       for ii=1:nPermIter
           xrand = allx(randperm(numel(allx)));
           yrand = ally(randperm(numel(ally)));
%            xrand = nan(size(allx));
%            yrand = nan(size(ally));
%            % shuffling within subject only
%            for ss=1:nSubj
%                inds = slabs==ss;
%                dat2shuff = allx(inds);
%                xrand(inds) = dat2shuff(randperm(numel(dat2shuff)));
%                dat2shuff = ally(inds);
%                yrand(inds) = dat2shuff(randperm(numel(dat2shuff)));
%            end
           rho = corr(xrand,yrand);
           shuff_corrs(vv,cc,ii) = rho;
       end
    end
end

% two-tailed p-value comparing real correlation to null correlations
p_corraccRT = 2*min(cat(3, mean(shuff_corrs>=repmat(real_corrs,1,1,nPermIter),3),  mean(shuff_corrs<=repmat(real_corrs,1,1,nPermIter),3)),[],3);

%%
if plotSpatAccVsRunRTBehav
    sub_colors = plasma(nSubj+1);
    biglims=[0,1];
    rtlims=[0,1.8];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            allx=[];ally=[];
            slabs=[];
               
            for ss=1:nSubj
                for run=1:nRuns
            
                    rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    
                    x = mean(sconf>0);
                    y = median(rts);
                    
                    allx=[allx;x];
                    ally=[ally;y];
                    slabs=[slabs;ss];
                    h=plot(x,y,'.','Color',sub_colors(ss,:));
                    
                end
    
%                 h=plot(means(:,1),means(:,2),'o','Color',sub_colors(ss,:));
                if run==1 
                    lh=[lh,h];
                end
                plot(mean(allx(slabs==ss)), mean(ally(slabs==ss)),'.','Color',sub_colors(ss,:),'MarkerSize',20)
            
            end
%             
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[0,1];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            
            rho = corr(allx,ally);
            my_p = p_corraccRT(vv,cc);
            set(gca, 'FontSize', 12)
            set(gca,'XLim',biglims,'YLim',rtlims)
%             plot([0.5, 0.5],rtlims,'-','Color','k');
%             plot([0.5,0.5],biglims,'-','Color','k');
            axis square
            if vv>9
                xlabel('Spatial decoding accuracy')
            end
            ylabel('RT (s)')
            

            title(sprintf('%s\nrho=%.2f, p=%.3f', ROI_names{spat_areas(vv)},rho,my_p));
        end

        suptitle(sprintf('Mean RT vs. spatial decoding in each vis area (runs) \n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end

%%
if plotSpatAccVsRunRTBehav
    sub_colors = plasma(nSubj+1);
    zlims=[-4,4];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            allx=[];ally=[];
            slabs=[];
               
            for ss=1:nSubj
                for run=1:nRuns
            
                    rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    
                    x = mean(sconf>0);
                    y = median(rts);
                    
                    allx=[allx;x];
                    ally=[ally;y];
                    slabs=[slabs;ss];
                   
                    
                end
                allx(slabs==ss) = zscore(allx(slabs==ss));
                ally(slabs==ss) = zscore(ally(slabs==ss));
                
                h=plot(allx(slabs==ss),ally(slabs==ss),'.','Color',sub_colors(ss,:));
%                 h=plot(means(:,1),means(:,2),'o','Color',sub_colors(ss,:));
                if run==1 
                    lh=[lh,h];
                end
                plot(mean(allx(slabs==ss)), mean(ally(slabs==ss)),'.','Color',sub_colors(ss,:),'MarkerSize',20)
            
            end
%             
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[0,1];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            
            rho = corr(allx,ally);
          
            set(gca, 'FontSize', 12)
            set(gca,'XLim',zlims,'YLim',zlims)
%             plot([0.5, 0.5],rtlims,'-','Color','k');
%             plot([0.5,0.5],biglims,'-','Color','k');
            axis square
            if vv>9
                xlabel('Spatial decoding accuracy')
            end
            ylabel('RT (s)')
            

            title(sprintf('%s\nrho=%.2f', ROI_names{spat_areas(vv)},rho));
        end

        suptitle(sprintf('Mean RT vs. spatial decoding in each vis area (runs) \n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end

%%
if plotSpatAccVsRunRTBoot
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;
            
            histvals = squeeze(shuff_corrs(vv,cc,:));
            histogram(histvals(:),100,'FaceColor',[0.8, 0.8, 0.8],'EdgeColor','none');
            xlim([-0.5, 0.5]);
            plot([real_corrs(vv,cc),real_corrs(vv,cc)], get(gca,'YLim'),'Color','k')
            my_p = p_corraccRT(vv,cc);
            title(sprintf('%s\np=%.3f', ROI_names{spat_areas(vv)},my_p));
        end

        suptitle(sprintf('Correlation between runwise RT and runwise decoding accuracy\n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end

%%
if plotSpatAccVsAllRTBehav
    rtlims=[0,1.8];
    biglims=[0,1];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            
            means=zeros(nSubj, 2);
            for ss=1:nSubj
                rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc));

                means(ss,1) = mean(sconf>0);
                means(ss,2) = median(rts);
            end

            h=plot(means(:,1),means(:,2),'o','Color',col(cc,:));
            rho=corr(means(:,1),means(:,2));
            % plot a best fit line
% %             [p,S]=polyfit(means(:,1),means(:,2),1);
% %             xeval=biglims;
% %             ypred=xeval*p(1)+p(2);
% %             plot(xeval,ypred,'-','Color','k')
            
            lh=[lh,h];
           

            set(gca, 'FontSize', 12)
            set(gca,'XLim',biglims,'YLim',rtlims)
            plot(biglims,[0,0],'-','Color','k');
            plot([0.5,0.5],rtlims,'-','Color','k');
            axis square
            if vv>9
                xlabel('spatial decoding')
            end
            ylabel('RT (s)')

            title(sprintf('%s\nrho=%.3f', ROI_names{spat_areas(vv)},rho));
        end

        suptitle(sprintf('Median RT vs spatial decoding each visual area\n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end

%%
if plotSpatAccVsAllAccBehav
    biglims=[0,1];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            
            means=zeros(nSubj, 2);
            for ss=1:nSubj
                correct=squeeze(behavcorrect_allsubs(ss,condlabs_allsubs(ss,:)==cc));
                sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc));

                means(ss,1) = mean(sconf>0);
                means(ss,2) = mean(correct);
            end

            h=plot(means(:,1),means(:,2),'o','Color',col(cc,:));
            
            rho = corr(means(:,1),means(:,2));
           
            lh=[lh,h];
           
            set(gca, 'FontSize', 12)
            set(gca,'XLim',biglims,'YLim',biglims)
            plot(biglims,[0.5, 0.5],'-','Color','k');
            plot([0.5,0.5],biglims,'-','Color','k');
            axis square
            if vv>9
                xlabel('spatial decoding')
            end
            ylabel('behav accuracy')

            title(sprintf('%s\nrho=%.3f', ROI_names{spat_areas(vv)},rho));
        end

        suptitle(sprintf('Behav accuracy vs spatial decoding each visual area\n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end

%%
if plotSpatAccVsRunAccBehav
    
    sub_colors = plasma(nSubj+1);
    xlims=[0,1];
    ylims=[0,1];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            allx=[];ally=[];
            slabs=[];
               
            for ss=1:nSubj
                for run=1:nRuns
            
                    correct=squeeze(behavcorrect_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run));
%                     rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    
                    x = mean(sconf>0);
                    y = mean(correct);
                    
                    allx=[allx;x];
                    ally=[ally;y];
                    slabs=[slabs;ss];
                    h=plot(x,y,'.','Color',sub_colors(ss,:));
                    
                end
    
%                 h=plot(means(:,1),means(:,2),'o','Color',sub_colors(ss,:));
                if run==1 
                    lh=[lh,h];
                end
                plot(mean(allx(slabs==ss)), mean(ally(slabs==ss)),'.','Color',sub_colors(ss,:),'MarkerSize',20)
            
            end
%             
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[0,1];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            
            rho = corr(allx,ally);
%             my_p = p_corraccRT(vv,cc);
            set(gca, 'FontSize', 12)
            set(gca,'XLim',xlims,'YLim',ylims)
%             plot([0.5, 0.5],rtlims,'-','Color','k');
%             plot([0.5,0.5],biglims,'-','Color','k');
            axis square
            if vv>9
                xlabel('Spatial decoding accuracy')
            end
            ylabel('Behavioral accuracy')
            

            title(sprintf('%s\nrho=%.2f', ROI_names{spat_areas(vv)},rho));
        end

        suptitle(sprintf('Mean behav acc vs. spatial decoding in each vis area (runs) \n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
    

end

%%

if plotSpatConfVsRT
    
    sub_colors = plasma(nSubj+1);
    rtlims=[0,3];
    conflims=[-40,40];
%     for ss=1:nSubj
%         for cc=[1]
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            allx = [];ally=[];
            allh=[];
            for ss=1:nSubj

                x=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));
                y=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));

                scatter(x,y,4,'MarkerFaceColor',[sub_colors(ss,:)],'MarkerFaceAlpha',0.2, 'MarkerEdgeColor','none');
                h=plot(mean(x),mean(y),'.','Color',sub_colors(ss,:),'MarkerSize',20);
                allh=[allh,h];
                
                allx=[allx;x];
                ally=[ally;y'];
                
            end
            
            for hi=1:length(allh)
                uistack(allh(hi),'top')                
            end
            
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[min(conflims),max(conflims)];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            rho = corr(allx,ally);
            set(gca, 'FontSize', 12)
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
%             set(gca,'XLim',conflims,'YLim',rtlims);
            axis square
            if vv>9                
                xlabel('Confidence of classifer (a.u.)')
            end
            ylabel('RT (sec)')
            
            title(sprintf('%s\nrho=%.2f', ROI_names{spat_areas(vv)},rho));
            
        end

        suptitle(sprintf('RT vs Classifier Confidence %s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
%     end
end

%%
if plotSpatConfVsRTZscore
    
    sub_colors = plasma(nSubj+1);
    rtlims=[0,3];
    conflims=[-40,40];
%     for ss=1:nSubj
%         for cc=[1]
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            allx = [];ally=[];
            allh=[];
            for ss=1:nSubj

                x=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));
                y=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));

                x=zscore(x);
                y=zscore(y);
                
                scatter(x,y,4,'MarkerFaceColor',[sub_colors(ss,:)],'MarkerFaceAlpha',0.2, 'MarkerEdgeColor','none');
                h=plot(mean(x),mean(y),'.','Color',sub_colors(ss,:),'MarkerSize',20);
                allh=[allh,h];
                
                allx=[allx;x];
                ally=[ally;y'];
                
            end
            
            for hi=1:length(allh)
                uistack(allh(hi),'top')                
            end
            
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[min(conflims),max(conflims)];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            rho = corr(allx,ally);
            set(gca, 'FontSize', 12)
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
%             set(gca,'XLim',conflims,'YLim',rtlims);
            axis square
            if vv>9                
                xlabel('Confidence of classifer (a.u.)')
            end
            ylabel('RT (sec)')
            
            title(sprintf('%s\nrho=%.2f', ROI_names{spat_areas(vv)},rho));
            
        end

        suptitle(sprintf('RT vs Classifier Confidence %s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
%     end
end
%%
sub_colors = plasma(nSubj+1);
rtlims=[0, 1.5];
conflims=[-25,25];
nConfBins=2;
if plotSpatConfVsRTBins
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            allx=[];ally=[];
            for ss=1:nSubj

                x=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));
                y=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1))';
                
                % binning the data into roughly equal sized bins using
                % percentiles
                bin_edges = prctile(x, linspace(0,100,nConfBins+1));        
                x_bins=nan(nConfBins,1);y_bins=nan(nConfBins,1);
                indsused=zeros(size(x));
                for bb=1:nConfBins
                    inds=(x>=bin_edges(bb) & x<bin_edges(bb+1)) | (bb==nConfBins & x==bin_edges(bb+1));
                    assert(~any(indsused(inds)==1))
                    indsused(inds) = 1;
                    x_bins(bb) = median(x(inds));
                    y_bins(bb) = median(y(inds));
                end
                assert(all(indsused==1))
              
                assert(all(x_bins<conflims(2) & x_bins>conflims(1)) && all(y_bins<rtlims(2) & y_bins>rtlims(1)))
                plot(x_bins,y_bins,'.','Color',sub_colors(ss,:));
                
                allx=[allx;x_bins];
                ally=[ally;y_bins];
            end
            set(gca, 'FontSize', 12)
            % plot a best fit line
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[min(conflims),max(conflims)];
%             ypred=xeval*p(1)+p(2);
%             plot(xeval,ypred,'-','Color','k')
%             sst = sum((ally-mean(ally)).^2);
%             sse = sum((ally-polyval(p,allx)).^2);
%             r2=1-(sse/sst);
            rho = corr(allx,ally);
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
            set(gca,'XLim',conflims,'YLim',rtlims);
            axis square
            if vv>9
                xlabel('Confidence of classifer (a.u.)')
            end
            ylabel('RT (sec)')
            
            title(sprintf('%s\nrho=%.2f', ROI_names{spat_areas(vv)},rho));
            
        end

        suptitle(sprintf('RT vs Classifier Confidence %s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
end

%%

if plotSpatConfVsCorrect
    conflims=[-5,5];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(spat_areas)));
        npy = ceil(length(spat_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;
            means = zeros(nSubj,2);
            for ss=1:nSubj

                x=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc));
                y=squeeze(behavcorrect_allsubs(ss,condlabs_allsubs(ss,:)==cc));

%                 plot(x,y,'.','Color',sub_colors(ss,:));
                means(ss,1) = mean(x(y==0));
                means(ss,2) = mean(x(y==1));
                plot(means(ss,:),[0,1], '-','Color',sub_colors(ss,:));
%                 plot(mean(x(y==1)),1, '.','Color',sub_colors(ss,:));
            end
            set(gca, 'FontSize', 12)
            plot(mean(means(:,1)),0,'o','Color','k');
            plot(mean(means(:,2)),1,'o','Color','k');
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
            set(gca,'XLim',conflims,'YLim',[-0.2, 1.2]);
            axis square
            if vv>9
                xlabel('Confidence of classifer (a.u.)')
            end
            ylabel('Correctness')
            
            title(sprintf('%s', ROI_names{spat_areas(vv)}));
            
        end

        suptitle(sprintf('Correctness vs Classifier Confidence\n%s condition',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
end