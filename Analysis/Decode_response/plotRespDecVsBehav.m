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

resp_areas = [12:14];  % visual ROIs
% resp_areas = [12:15,22];

nROIs_resp = length(resp_areas);
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

plotRespAccOverTime = 1;

plotRespAccVsAllAccBehav=1;
plotRespAccVsRunAccBehav=1;

plotRespAccVsAllRTBehav=1;
plotRespAccVsRunRTBehav=1;
plotRespAccVsRunRTBehavZscore=1;
plotRespAccVsRunRTBoot=1;

plotRespConfVsRT = 1;
plotRespConfVsRTZscore=1;
plotRespConfVsRTBins = 1;
plotRespConfVsCorrect = 1;
%% load results
nTrialsTotal = 400;
nRuns=20;
nTrialsPerRun=20;
nTRs_out = 30;
trDur = 0.8;
tax = trDur*(0:nTRs_out-1);
evts2plot = [3.5, 4.5, 16.5, 18.5];

respacc_allsubs = nan(nSubj,nROIs_resp,nConds);
respacc_overtime_allsubs = nan(nSubj,nROIs_resp,nConds, nTRs_out);
% respacc_allsubs = nan(nSubj,nROIs_resp,nConds);
respconf = nan(nSubj, nROIs_resp, nTrialsTotal);
respconf_overtime = nan(nSubj, nROIs_resp, nTrialsTotal, nTRs_out);
% respconf = nan(nSubj, nROIs_resp, nTrialsTotal);
rt_allsubs = nan(nSubj, nTrialsTotal);
behavcorrect_allsubs = nan(nSubj, nTrialsTotal);
condlabs_allsubs = nan(nSubj, nTrialsTotal);

runlabs=repelem(1:nRuns, nTrialsPerRun);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    % loading performance of resp position classifier
    save_dir = fullfile(exp_path,'Analysis','Decode_response','Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    
   
    
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    fn2load = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load,'rt','correct','condlabs');
   
    
    respacc_allsubs(ss,:,:) = squeeze(allacc(resp_areas,:));
    respconf(ss,:,:) = allconf(resp_areas,:);
    rt_allsubs(ss,:) = rt;
    behavcorrect_allsubs(ss,:) = correct;
    condlabs_allsubs(ss,:) = condlabs;
    clear allacc
    clear allconf
    
    save_dir = fullfile(exp_path,'Analysis','Decode_response','Decoding_results');
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    
    respacc_overtime_allsubs(ss,:,:,:) = squeeze(allacc(resp_areas,:,:));
    respconf_overtime(ss,:,:,:,:) = squeeze(allconf(resp_areas,:,:,:));
    
    % loading performance of resp classifier
%     save_dir = fullfile(exp_path,'Analysis','Decode_resp','Decoding_results');
%     fn2load = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
%     load(fn2load);
%     
%     respacc_allsubs(ss,:,:) = allacc(resp_areas,:);
%     respconf(ss,:,:) = allconf(resp_areas,:);
%     clear allacc
%     clear allconf
end

assert(~any(isnan(respacc_allsubs(:))))
% assert(~any(isnan(respacc_allsubs(:))))


%% Plot timecourse of decoding acc, for different bins of RTs
% combining all subjects - concatenate all trials from all subjects, then
% break into roughly equal sized bins.
nRTbins = 3;

rt_colors = parula(nRTbins+1);
lw=1;
fs = 10;
if plotRespAccOverTime
    
   
    for cc = 1:nConds
        nplots_motor = ceil(sqrt(nROIs_resp));
        figure();hold all;

        for vi = 1:nROIs_resp

            subplot(nplots_motor,ceil(nROIs_resp/nplots_motor),vi);hold all;

            lh=[];

            rtvals=[];
            confvals = [];
            for ss = 1:nSubj

                rtvals = [rtvals;squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1))'];
                confvals = [confvals;squeeze(respconf_overtime(ss,vi,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1,:))];    
                
            end
            % binning the data into roughly equal sized bins using
            % percentiles
            bin_edges = prctile(rtvals, linspace(0,100,nRTbins+1));        
  
            dec_bins=nan(nRTbins,nTRs_out);
            confmean_bins = nan(nRTbins, nTRs_out);
            confsem_bins = nan(nRTbins, nTRs_out);
            
            indsused=zeros(size(rtvals));
            rtlabs = [];
            for bb=1:nRTbins
                inds=(rtvals>=bin_edges(bb) & rtvals<bin_edges(bb+1)) | (bb==nRTbins & rtvals==bin_edges(bb+1));
                assert(~any(indsused(inds)==1))
                indsused(inds) = 1;
                
                % converting from confidence to acc here
                dec_bins(bb,:) = mean(confvals(inds,:)>0,1);
                confmean_bins(bb,:) = mean(confvals(inds,:),1);
                confsem_bins(bb,:) = std(confvals(inds,:),1)./sqrt(sum(inds)-1);
                
                rtlabs{bb} = sprintf('RT %.2f - %.2f sec',min(rtvals(inds)),max(rtvals(inds)));
            end
            assert(all(indsused==1))

        
           
            for bb = 1:nRTbins
                
                meanvals = dec_bins(bb,:)';
                semvals = nan(size(meanvals));
%                 vals = squeeze(dec_bins(:,bb,:));
%                 if nSubj==1
%                     meanvals =squeeze(vals)';
%                     semvals = [];
%                 else
%                     meanvals = squeeze(mean(vals,1))';
%                     semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
%                 end

                lh=[lh,plot(tax,meanvals,'-','Color',rt_colors(bb,:),'LineWidth',lw)];

                bandedError_MMH(tax, meanvals',semvals',rt_colors(bb,:), 0.5);
            end
%             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
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

            if vi==nROIs_resp
                legend(lh,rtlabs,'FontSize', fs);
            end

            title(sprintf('%s', ROI_names{resp_areas(vi)}));
        end
        suptitle(sprintf('%s',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
              
    end
end

%% Plot timecourse of decoding confidence, for different bins of RTs
% combining all subjects - concatenate all trials from all subjects, then
% break into roughly equal sized bins.
nRTbins = 3;
conflims = [-1,5];
rt_colors = parula(nRTbins+1);
lw=1;
fs = 10;
if plotRespAccOverTime
    
   
    for cc = 1:nConds
        nplots_motor = ceil(sqrt(nROIs_resp));
        figure();hold all;

        for vi = 1:nROIs_resp

            subplot(nplots_motor,ceil(nROIs_resp/nplots_motor),vi);hold all;

            lh=[];

            rtvals=[];
            confvals = [];
            for ss = 1:nSubj

                rtvals = [rtvals;squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1))'];
                confvals = [confvals;squeeze(respconf_overtime(ss,vi,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1,:))];    
                
            end
            % binning the data into roughly equal sized bins using
            % percentiles
            bin_edges = prctile(rtvals, linspace(0,100,nRTbins+1));        
  
            dec_bins=nan(nRTbins,nTRs_out);
            confmean_bins = nan(nRTbins, nTRs_out);
            confsem_bins = nan(nRTbins, nTRs_out);
            
            indsused=zeros(size(rtvals));
            rtlabs = [];
            for bb=1:nRTbins
                inds=(rtvals>=bin_edges(bb) & rtvals<bin_edges(bb+1)) | (bb==nRTbins & rtvals==bin_edges(bb+1));
                assert(~any(indsused(inds)==1))
                indsused(inds) = 1;
                
                % converting from confidence to acc here
                dec_bins(bb,:) = mean(confvals(inds,:)>0,1);
                confmean_bins(bb,:) = mean(confvals(inds,:),1);
                confsem_bins(bb,:) = std(confvals(inds,:),1)./sqrt(sum(inds)-1);
                
                rtlabs{bb} = sprintf('RT %.2f - %.2f sec',min(rtvals(inds)),max(rtvals(inds)));
            end
            assert(all(indsused==1))

        
           
            for bb = 1:nRTbins
                
                meanvals = confmean_bins(bb,:)';
                semvals = confsem_bins(bb,:)';
%                 vals = squeeze(dec_bins(:,bb,:));
%                 if nSubj==1
%                     meanvals =squeeze(vals)';
%                     semvals = [];
%                 else
%                     meanvals = squeeze(mean(vals,1))';
%                     semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
%                 end

                lh=[lh,plot(tax,meanvals,'-','Color',rt_colors(bb,:),'LineWidth',lw)];

                bandedError_MMH(tax, meanvals',semvals',rt_colors(bb,:), 0.5);
            end
%             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
            set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
            ylim(conflims);
            plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
            for ee = 1:length(evts2plot)
                plot([evts2plot(ee),evts2plot(ee)],conflims,'-','Color',[0.8, 0.8, 0.8]);
            end

            if vi==1
                xlabel('Time(s)');
                ylabel('Classifier "confidence"');
            end

            if vi==nROIs_resp
                legend(lh,rtlabs,'FontSize', fs);
            end

            title(sprintf('%s', ROI_names{resp_areas(vi)}));
        end
        suptitle(sprintf('%s',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1200,800]);
              
    end
end

%% Plot timecourse of decoding, for different bins of RTs
% Single subjects, very noisy
nRTbins = 3;
rtlabs = {'rt bin 1 (shortest)','rt bin 2','rt bin 3 (longest)'};
nperbin = zeros(nRTbins, nSubj, nConds);
bin_edges_all = zeros(nSubj, nConds, nRTbins+1);
rt_colors = parula(nRTbins+1);
lw=1;
acclims = [0.2, 1.0];
fs = 10;
if plotRespAccOverTime
    
    
    for cc = 1:nConds
%             nplots_motor = ceil(sqrt(nROIs_resp));
        figure();hold all;
        pi=0;
        for ss = 1:nSubj
            for vi = 1:nROIs_resp
                pi=pi+1;
                subplot(nSubj,nROIs_resp,pi);hold all;

                lh=[];

                rtvals = squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1))';
                confvals = squeeze(respconf_overtime(ss,vi,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1,:));    
                % binning the data into roughly equal sized bins using
                % percentiles
                bin_edges = prctile(rtvals, linspace(0,100,nRTbins+1));
                bin_edges_all(ss,cc,:) = bin_edges;
                rt_bins=nan(nRTbins,1);dec_bins=nan(nRTbins,nTRs_out);
                indsused=zeros(size(rtvals));
%                 rtlabs = [];
                for bb=1:nRTbins
                    inds=(rtvals>=bin_edges(bb) & rtvals<bin_edges(bb+1)) | (bb==nRTbins & rtvals==bin_edges(bb+1));
                    assert(~any(indsused(inds)==1))
                    indsused(inds) = 1;
                    rt_bins(bb) = median(rtvals(inds));
                    nperbin(bb,ss,cc) = sum(inds);
                    % converting from confidence to acc here
                    dec_bins(bb,:) = mean(confvals(inds,:)>0,1);
%                     rtlabs{bb} = sprintf('RT %.2f - %.2f sec',min(rtvals(inds)),max(rtvals(inds)));
                end
                assert(all(indsused==1))
                
                
                for bb = 1:nRTbins

                    lh=[lh,plot(tax,dec_bins(bb,:),'-','Color',rt_colors(bb,:),'LineWidth',lw)];

%                     bandedError_MMH(tax, dec_bins(bb,:),[],squeeze(col_conds(cc,:,bb)), 0.5);
                end
    %             errorbar(tax, meanvals, semvals,'Color',col(cc,:),'LineWidth',lw);
                set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
                ylim(acclims);
                plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
                for ee = 1:length(evts2plot)
                    plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
                end

                if vi==1
                    xlabel('Time(s)');
                    ylabel(sprintf('S%02d',sublist(ss)));
                end

                if vi==nROIs_resp && ss==1
                    legend(lh,rtlabs,'FontSize', fs);
                end

                if ss==1
                    title(sprintf('%s', ROI_names{resp_areas(vi)}));
                end
            end
           
        end
        suptitle(sprintf('%s',condLabStrs{cc}));
        set(gcf,'Color','w')
        set(gcf,'Position',[200,200,1600,1200]);

    end
end


%% checking to make sure the confidence values give the correct accuracy calculations
% note they aren't perfectly identical (but should be very close) due to
% the way that I averaged over groups to get accuracy and then averaged
% those accuracy values to get teh final accuracy (what is stored in allacc
% variable) 
acc_check = nan(nSubj, numel(resp_areas), nConds);
for ss=1:nSubj
    for vv=1:numel(resp_areas)
        for cc=1:nConds
            conf = respconf(ss,vv,condlabs_allsubs(ss,:)==cc);
            acc_check(ss,vv,cc) = mean(conf>0);
        end
    end
end
vals = acc_check;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));


%% Bin the "confidence" values, and plot versus RT
sub_colors = plasma(nSubj+1);
rtlims=[0, 1.5];
conflims=[-40,40];
nConfBins=4;
if plotRespConfVsRTBins
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(resp_areas)));
        npy = ceil(length(resp_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(resp_areas)
            subplot(npx,npy,vv);hold all;

            allx=[];ally=[];
            for ss=1:nSubj

                x=squeeze(respconf(ss,vv,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));
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
            if vv>2
                xlabel('Confidence of classifer (a.u.)')
            end
            ylabel('RT (sec)')
            
            title(sprintf('%s\nrho=%.2f', ROI_names{resp_areas(vv)},rho));
            
        end

        suptitle(sprintf('RT vs Classifier Confidence %s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
end



%% bootstrap correlations of RT versus spat decoding performance (using run-averaged values)
nPermIter=1000;
shuff_corrs = nan(nROIs_resp, nConds, nPermIter);
real_corrs = nan(nROIs_resp, nConds);

for vv = 1:nROIs_resp
    for cc=1:nConds
        % for this ROI and cond, want to get a single value for each
        % subject and each run of each task (e.g. 20*6 total)
        allx=[];ally=[];
        slabs=[];
        for run=1:nRuns
            means=zeros(nSubj, 2);
            for ss=1:nSubj
                rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                sconf=squeeze(respconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                
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
if plotRespAccVsRunRTBehav
    sub_colors = plasma(nSubj+1);
    biglims=[0,1];
    rtlims=[0,1.8];
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(resp_areas)));
        npy = ceil(length(resp_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(resp_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            
            allx=[];ally=[];
            slabs=[];
               
            for ss=1:nSubj
                for run=1:nRuns
            
                    rts=squeeze(rt_allsubs(ss,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    sconf=squeeze(respconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run & behavcorrect_allsubs(ss,:)==1));
                    
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
            if vv>2
                xlabel('resp decoding accuracy')
            end
            ylabel('RT (s)')
            

            title(sprintf('%s\nrho=%.2f, p=%.3f', ROI_names{resp_areas(vv)},rho,my_p));
        end

        suptitle(sprintf('Mean RT vs. resp decoding in each vis area (runs) \n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end


%%
if plotRespAccVsRunRTBoot
    for cc=1:nConds
        figure;hold all;
        npx = ceil(sqrt(length(resp_areas)));
        npy = ceil(length(resp_areas)/npx);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(resp_areas)
            subplot(npx,npy,vv);hold all;
            
            histvals = squeeze(shuff_corrs(vv,cc,:));
            histogram(histvals(:),100,'FaceColor',[0.8, 0.8, 0.8],'EdgeColor','none');
            xlim([-0.5, 0.5]);
            plot([real_corrs(vv,cc),real_corrs(vv,cc)], get(gca,'YLim'),'Color','k')
            my_p = p_corraccRT(vv,cc);
            title(sprintf('%s\np=%.3f', ROI_names{resp_areas(vv)},my_p));
        end

        suptitle(sprintf('Correlation between runwise RT and runwise decoding accuracy\n%s',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end
