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
resp_areas = [12:14];

nROIs_spat = length(spat_areas);
nROIs_resp = length(resp_areas);

nVox2Use = 10000;
nPermIter=1000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

acclims = [0.4, 0.9];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);


condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);
plotAccOverall=1;
plotAccEachRun=1;
plotConfTrials=1;
plotConfBins=1;

nConfBins=8;

%% load results
nTrialsTotal = 400;
nRuns=20;
nTrialsPerRun=20;

spatacc_allsubs = nan(nSubj,nROIs_spat,nConds);
respacc_allsubs = nan(nSubj,nROIs_resp,nConds);
spatconf = nan(nSubj, nROIs_spat, nTrialsTotal);
respconf = nan(nSubj, nROIs_resp, nTrialsTotal);
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
    assert(size(allacc,1)==numel(ROI_names));
    
    spatacc_allsubs(ss,:,:) = mean(squeeze(allacc(spat_areas,:,:)),3);
    spatconf(ss,:,:) = allconf(spat_areas,:);
    rt_allsubs(ss,:) = rt;
    behavcorrect_allsubs(ss,:) = correct;
    condlabs_allsubs(ss,:) = condlabs;
    clear allacc
    clear allconf
    
    % loading performance of response classifier
    save_dir = fullfile(exp_path,'Analysis','Decode_response','Decoding_results');
%     fn2load = fullfile(save_dir,sprintf('ClassifyActualResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    respacc_allsubs(ss,:,:) = allacc(resp_areas,:);
    respconf(ss,:,:) = allconf(resp_areas,:);
    clear allacc
    clear allconf
end

assert(~any(isnan(spatacc_allsubs(:))))
assert(~any(isnan(respacc_allsubs(:))))

%% checking to make sure the confidence values give the correct accuracy calculations
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

% plot_barsAndStars(meanvals,semvals,[],...
%         [],0.5,acclims,ROI_names(spat_areas),condLabStrs,...
%         'Accuracy','Spatial Memory Position',col)


%%

if plotAccOverall
    
    for rr = 1:3
         
        figure;hold all;
        npy = ceil(sqrt(length(spat_areas)));
        npx = ceil(length(spat_areas)/npy);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            for cc=1:nConds
                x=squeeze(respacc_allsubs(:,rr,cc));
                y=squeeze(spatacc_allsubs(:,vv,cc));

                lh=[lh,plot(x,y,'o','Color',col(cc,:))];
                
                linefit = cat(1,[x, ones(nSubj,1)])\y;
                yhat = x * linefit(1) + linefit(2);
                sse = sum((yhat-y).^2);
                sst = sum((y-mean(y)).^2);
                r2 = 1-sse/sst;
                rho = corr(x,y);
                xpts = [min(x), max(x)];
                ypts = xpts * linefit(1) + linefit(2);
                plot(xpts, ypts, '-', 'Color',col(cc,:))
                
            end

            set(gca, 'FontSize', 12)
            set(gca,'XLim',acclims,'YLim',acclims)
            axis square
            if vv==9
                xlabel('response decoding')
                ylabel('spatial decoding')
            else
                xticks([])
                yticks([])
            end
            title(sprintf('%s', ROI_names{spat_areas(vv)}));
            if vv==numel(spat_areas)
                legend(lh,condLabStrs)
            end
        end

        suptitle(sprintf('Spat decoding in each visual area vs. response decoding in %s',ROI_names{resp_areas(rr)}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')

    end
    
%     rr=1;
%     ROI_names{resp_areas(rr)} = 'S1/M1/PMc avg';

    figure;hold all;
    npy = ceil(sqrt(length(spat_areas)));
    npx = ceil(length(spat_areas)/npy);
    % plot each ROI as a separate subplot.

    for vv = 1:numel(spat_areas)
        subplot(npx,npy,vv);hold all;

        lh=[];
        for cc=1:nConds
            x=squeeze(mean(respacc_allsubs(:,:,cc),2));
            y=squeeze(spatacc_allsubs(:,vv,cc));
            
            lh=[lh,plot(x,y,'o','Color',col(cc,:))];
            
            linefit = cat(1,[x, ones(nSubj,1)])\y;
            yhat = x * linefit(1) + linefit(2);
            sse = sum((yhat-y).^2);
            sst = sum((y-mean(y)).^2);
            r2 = 1-sse/sst;
            rho = corr(x,y);
            xpts = [min(x), max(x)];
            ypts = xpts * linefit(1) + linefit(2);
            plot(xpts, ypts, '-', 'Color',col(cc,:))
        end
        
        set(gca, 'FontSize', 12)
        set(gca,'XLim',acclims,'YLim',acclims)
        axis square
        if vv==9
            xlabel('response decoding')
            ylabel('spatial decoding')
        else
            xticks([])
            yticks([])
        end
        title(sprintf('%s', ROI_names{spat_areas(vv)}));
        if vv==numel(spat_areas)
            legend(lh,condLabStrs)
        end
    end

    suptitle(sprintf('Spat decoding in each visual area vs. response decoding in %s','average M1/S1/PMC'))
    set(gcf,'Position',[200,200,1400,1400]);

    set(gcf,'Color','w')

end

%%
if plotAccEachRun
    biglims=[0,1];
    for cc=1:nConds
        figure;hold all;
        npy = ceil(sqrt(length(spat_areas)));
        npx = ceil(length(spat_areas)/npy);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            lh=[];
            allvals = [];
            for run=1:nRuns
                means=zeros(nSubj, 2);
                for ss=1:nSubj
                    rconf=mean(squeeze(respconf(ss,:,condlabs_allsubs(ss,:)==cc & runlabs==run)), 1);
                    sconf=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & runlabs==run));
                    means(ss,1) = mean(rconf>0);
                    means(ss,2) = mean(sconf>0);
                end
                allvals = [allvals; means];
                h=plot(means(:,1),means(:,2),'o','Color',col(cc,:));
               
                if run==1 
                    lh=[lh,h];
                end
            end
            
            x = allvals(:,1);
            y = allvals(:,2);

            linefit = cat(1,[x, ones(size(x,1),1)])\y;
            yhat = x * linefit(1) + linefit(2);
            sse = sum((yhat-y).^2);
            sst = sum((y-mean(y)).^2);
            r2 = 1-sse/sst;
            [rho,p] = corr(x,y);
            xpts = [min(x), max(x)];
            ypts = xpts * linefit(1) + linefit(2);
            plot(xpts, ypts, '-', 'Color','k')

            
%             zx = zscore(x);
%             zy = zscore(y);
%             zlinefit = cat(1,[zx, ones(size(zx,1),1)])\zy;
%             
%           

            set(gca, 'FontSize', 12)
            set(gca,'XLim',biglims,'YLim',biglims)
            plot(biglims,[0.5, 0.5],'-','Color','k');
            plot([0.5,0.5],biglims,'-','Color','k');
            axis square
            
            if vv==9
                xlabel('response decoding')
                ylabel('spatial decoding')
            else
                xticks([])
                yticks([])
            end

            title(sprintf('%s\nr=%.2f, p=%.2f', ROI_names{spat_areas(vv)}, rho, p));
        end

        suptitle(sprintf('Spat decoding in each visual area vs. response decoding in %s\n%s','avg M1/S1/PMC',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end

end


%%
rr=1;
sub_colors = plasma(nSubj+1);
conflims=[-40,40];
if plotConfTrials
% for ss=1:nSubj
    for cc=1:nConds
        figure;hold all;
        npy = ceil(sqrt(numel(spat_areas)));
        npx = ceil(numel(spat_areas)/npy);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

            allx=[];ally=[];
            for ss=1:nSubj

                x=squeeze(mean(respconf(ss,:,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1), 2));
                y=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc & behavcorrect_allsubs(ss,:)==1));
                
                plot(x,y,'.','Color',sub_colors(ss,:));
                
                allx=[allx;x];
                ally=[ally;y];
            end
           
            x = allx;
            y = ally;
            linefit = cat(1,[x, ones(size(x,1),1)])\y;
            yhat = x * linefit(1) + linefit(2);
            sse = sum((yhat-y).^2);
            sst = sum((y-mean(y)).^2);
            r2 = 1-sse/sst;
            [rho,p] = corr(x,y);
            xpts = [min(x), max(x)];
            ypts = xpts * linefit(1) + linefit(2);
            plot(xpts, ypts, '-', 'Color','k')

            set(gca, 'FontSize', 12)
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
            set(gca,'XLim',conflims,'YLim',conflims);
            plot(conflims,[0,0],'-','Color','k');
            plot([0,0],conflims,'-','Color','k');
            axis square
%             xticks([])
%             yticks([])
            if vv==9
                xlabel('Resp conf.')
                ylabel('Spatial conf.')
            else
                xticks([])
                yticks([])
            end
            
            title(sprintf('%s\nrho=%.2f, p=%.2f', ROI_names{spat_areas(vv)},rho,p));
            
        end

        suptitle(sprintf('Spat decoding in each visual area vs. response decoding in %s\n%s','avg M1/S1/PMc',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
% end
end

%%
rr=1;
sub_colors = plasma(nSubj+1);
conflims=[-10,30];
nConfBins=8;

if plotConfBins
    for cc=1:nConds
        figure;hold all;
        npy = ceil(sqrt(numel(spat_areas)));
        npx = ceil(numel(spat_areas)/npy);
        % plot each ROI as a separate subplot.

        for vv = 1:numel(spat_areas)
            subplot(npx,npy,vv);hold all;

%             for ss=1
            allx=[];ally=[];
            for ss=1:nSubj

                x=squeeze(mean(respconf(ss,:,condlabs_allsubs(ss,:)==cc),2));
                y=squeeze(spatconf(ss,vv,condlabs_allsubs(ss,:)==cc));
                
                [ysort,sorder] = sort(y, 'ascend');
                xbinned = reshape(x(sorder),nTrialsTotal/2/nConfBins, nConfBins);
                ybinned = reshape(ysort,nTrialsTotal/2/nConfBins, nConfBins);
                x_bins = nanmedian(xbinned,1);
                y_bins = nanmedian(ybinned,1);
               
                plot(x_bins,y_bins,'.','Color',sub_colors(ss,:));
                assert(min(x_bins)>conflims(1) & max(x_bins)<conflims(2))
                
                allx=[allx;x_bins'];
                ally=[ally;y_bins'];
            end
           
            x = allx;y = ally;
            linefit = cat(1,[x, ones(size(x,1),1)])\y;
            yhat = x * linefit(1) + linefit(2);
            sse = sum((yhat-y).^2);
            sst = sum((y-mean(y)).^2);
            r2 = 1-sse/sst;
            [rho,p] = corr(x,y);
            xpts = [min(x), max(x)];
            ypts = xpts * linefit(1) + linefit(2);
            plot(xpts, ypts, '-', 'Color','k')
%             rho = corr(allx,ally);
            % plot a best fit line
%             [p,S]=polyfit(allx,ally,1);
%             xeval=[min(conflims),max(conflims)];
%             pred=xeval*p(1)+p(2);
%             plot(xeval,pred,'-','Color','k')
            
            set(gca, 'FontSize', 12)
%             set(gca,'XTick',[0:90:360],'YTick',[0:90:360]);
%             set(gca,'XLim',conflims,'YLim',conflims);
    
%             plot(conflims,[0,0],'-','Color','k');
%             plot([0,0],conflims,'-','Color','k');
            plot(get(gca,'XLim'),[0,0],'-','Color','k');
            plot([0,0],get(gca,'YLim'),'-','Color','k');
            axis square
            xticks([])
            yticks([])
            if vv==9
                xlabel('Resp conf.')
                ylabel('Spatial conf.')
            else
                xticks([])
                yticks([])
            end
            title(sprintf('%s\nrho=%.2f, p=%.2f', ROI_names{spat_areas(vv)},rho,p));
            
        end

        suptitle(sprintf('Spat decoding in each visual area vs. response decoding in %s (binned)\n%s','avg M1/S1/PMc',condLabStrs{cc}))
        set(gcf,'Position',[200,200,1400,1400]);

        set(gcf,'Color','w')
    end
end 
