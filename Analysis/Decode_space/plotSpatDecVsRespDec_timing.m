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
     'S1','M1','Premotor',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

spat_areas = [1:5,10,11,6:9];  % visual ROIs
% spat_areas=[1:5];
resp_areas = [12:14];

nROIs_spat = length(spat_areas);
nROIs_resp = length(resp_areas);

nVox2Use = 10000;

% class_str = 'svmtrain_lin';
class_str = 'normEucDist';

lw=2;
fs=12;
ms=10;
acclims = [0.4, 0.9];
accdifflims = [-0.2, 0.2];
dprimelims = [-0.2, 1.4];
col = plasma(5);
col = col(2:2:end-1,:);

chance_val=0.5;
trDur = 0.8;
nTRs_out=30;
tax = trDur*(0:nTRs_out-1);
evts2plot = [3.5, 4.5, 16.5, 18.5];
% what is the first TR to look for a divergence between random and
% predictable conditions (this is the preview disk onset time)
min_ind = find(tax>evts2plot(1),1);

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);
plotSpatAccDiff = 1;
plotRespAcc = 1;
plotSpatDiffSS=1;

nConfBins=8;
nTrialsTotal = 400;
nRuns=20;
nTrialsPerRun=20;

nPermIter=1000;

%% load results


spatacc_allsubs = nan(nSubj,nROIs_spat,nConds, nTRs_out);
respacc_allsubs = nan(nSubj,nROIs_resp,nConds, nTRs_out);
respaccrand_allsubs = nan(nSubj,nROIs_resp,nConds, nTRs_out, nPermIter);
% spatconf = nan(nSubj, nROIs_spat, nTrialsTotal);
% respconf = nan(nSubj, nROIs_resp, nTrialsTotal);
rt_allsubs = nan(nSubj, nTrialsTotal);
behavcorrect_allsubs = nan(nSubj, nTrialsTotal);
condlabs_allsubs = nan(nSubj, nTrialsTotal);

runlabs=repelem(1:nRuns, nTrialsPerRun);

for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    % loading performance of spatial position classifier
    save_dir = fullfile(exp_path,'Analysis','Decode_space','Decoding_results');
    fn2load = fullfile(save_dir,sprintf('TrnSWMLoc_TestWM_TRbyTR_%s_max%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    
    spatacc_allsubs(ss,:,:,:) = mean(squeeze(allacc(spat_areas,:,:,:)),3);
%     spatconf(ss,:,:) = allconf(spat_areas,:);
    rt_allsubs(ss,:) = rt;
    behavcorrect_allsubs(ss,:) = correct;
    condlabs_allsubs(ss,:) = condlabs;
    clear allacc
    clear allconf
    
    % loading performance of response classifier
    save_dir = fullfile(exp_path,'Analysis','Decode_response','Decoding_results');
%     fn2load = fullfile(save_dir,sprintf('ClassifyActualResponse_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    fn2load = fullfile(save_dir,sprintf('ClassifyResponse_TRbyTR_%s_%dvox_%s.mat',class_str,nVox2Use,substr));
    load(fn2load);
    assert(size(allacc,1)==numel(ROI_names));
    respacc_allsubs(ss,:,:,:) = allacc(resp_areas,:,:);
    respaccrand_allsubs(ss,:,:,:,:) = allacc_rand(resp_areas,:,:,:);
%     respconf(ss,:,:) = allconf(resp_areas,:);
    clear allacc
    clear allconf
end

assert(~any(isnan(spatacc_allsubs(:))))
assert(~any(isnan(respacc_allsubs(:))))

%% now doing pairwise condition comparisons - sign-rank test

p_diff_spat = nan(numel(spat_areas),nTRs_out);
dir_diff = nan(numel(spat_areas),nTRs_out);
spat_diff_onsets = nan(numel(spat_areas),1);
for vv=1:numel(spat_areas)
    for tr=1:nTRs_out
        realvals = squeeze(spatacc_allsubs(:,vv,:,tr));
        [p,h,stats] = signrank(realvals(:,1),realvals(:,2));
        dir_diff(vv,tr) = double(mean(realvals(:,1))<mean(realvals(:,2)));
%         [h,p,ci,stats] = ttest(realvals(:,1),realvals(:,2));
        p_diff_spat(vv,tr) = p;
    end
    is_sig = p_diff_spat(vv,:)<0.05;
    is_sig(1:min_ind-1) = 0;
    % find the time of the TR at which the conditions diverged.
    if sum(is_sig)>0
        spat_diff_onsets(vv) = tax(find(is_sig & dir_diff(vv,:)==1,1));
    end
end
spat_diff_is_sig = p_diff_spat<0.05;    

%% compare response decoding against chance
% for each permutation iteration, use this test to compare real data for all subj to
% shuffled data for all subj.
resp_iters_sr = nan(numel(resp_areas), nConds, nTRs_out, nPermIter); 

for vv=1:numel(resp_areas)
    for cc=1:nConds
        for tr = 1:nTRs_out
            x = respacc_allsubs(:,vv,cc,tr);
        
            for ii=1:nPermIter
                y = respaccrand_allsubs(:,vv,cc,tr,ii);
                % compare the real values against the null, for this iteration.
                % w>0 means real>null, w<0 means real<null, w=0 means equal
                resp_iters_sr(vv,cc,tr,ii) = signrank_MMH(x,y);

            end
        end
    end
end

% final p value is the proportion of iterations where null was at least as
% large as the real (e.g. the test stat was 0 or negative)
p_resp = mean(resp_iters_sr<=0, 4);
resp_is_sig = p_resp<0.05;    % one tailed test

resp_sig_onsets = nan(numel(resp_areas),1);
for vv=1:numel(resp_areas)
    is_sig = p_resp(vv,1,:)<0.05;
    is_sig(1:min_ind-1) = 0;
    if sum(is_sig)>0
        resp_sig_onsets(vv) = tax(find(is_sig,1)); 
    end
end

%%
if plotSpatAccDiff
   
    sig_heights = linspace(0.15, 0.20,numel(spat_areas));

    area_colors = plasma(numel(spat_areas)+1);
    figure();hold all;
    
    lh=[];
    for vi = 1:numel(spat_areas)
 
        
       
        vals = squeeze(spatacc_allsubs(:,vi,:,:));
        % take random - predictable difference 
        vals = squeeze(vals(:,2,:)-vals(:,1,:));
        if nSubj==1
            meanvals =squeeze(vals)';
            semvals = [];
        else
            meanvals = squeeze(mean(vals,1))';
            semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
        end
        lh=[lh,plot(tax,meanvals,'-','Color',area_colors(vi,:),'LineWidth',lw)];
%         bandedError_MMH(tax, meanvals',semvals', area_colors(vi,:), 0.5);
       
        plot(spat_diff_onsets(vi), sig_heights(vi),'.','Color',area_colors(vi,:),'MarkerSize',ms);
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
%         ylim(acclims);
        plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],accdifflims,'-','Color',[0.8, 0.8, 0.8]);
        end
                
        xlabel('Time(s)');
        ylabel('Accuracy');
        
        title(sprintf('Random - Predictable'));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
    legend(lh,ROI_names(spat_areas));
end

%%
sig_heights = linspace(0.8, 0.82,numel(resp_areas));
if plotRespAcc
   
    area_colors = plasma(numel(resp_areas)+1);
    figure();hold all;
    
    lh=[];
    for vi = 1:numel(resp_areas)
 
        
        % predictable cond only
        vals = squeeze(respacc_allsubs(:,vi,1,:));
       
        if nSubj==1
            meanvals =squeeze(vals)';
            semvals = [];
        else
            meanvals = squeeze(mean(vals,1))';
            semvals = squeeze(std(vals,[],1))'./sqrt(nSubj);
        end
        lh=[lh,plot(tax,meanvals,'-','Color',area_colors(vi,:),'LineWidth',lw)];
%         bandedError_MMH(tax, meanvals',semvals', area_colors(vi,:), 0.5);
       
        plot(resp_sig_onsets(vi), sig_heights(vi),'.','Color',area_colors(vi,:),'MarkerSize',ms);
        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(acclims);
        plot(get(gca,'XLim'),[chance_val,chance_val],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],acclims,'-','Color',[0.8, 0.8, 0.8]);
        end
                
        xlabel('Time(s)');
        ylabel('Accuracy');
        
        title(sprintf('Predictable cond response decoding'));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1200,800]);
    legend(lh,ROI_names(resp_areas));
end


%%
if plotSpatDiffSS
    sub_colors=viridis(nSubj+1);
    nplots_vis = ceil(sqrt(numel(spat_areas)));
    figure();hold all;
    lims=[-0.5, 0.5];
    for vi = 1:numel(spat_areas)
       
        subplot(nplots_vis,ceil(numel(spat_areas)/nplots_vis),vi);hold all;

        cc=1;
        valsplot = spatacc_allsubs(:,vi,2,:) - spatacc_allsubs(:,vi,1,:);
        lh=[];
        substrs=[];
        for ss=1:nSubj
            lh=[lh, plot(tax,squeeze(valsplot(ss,:)),'-','Color',sub_colors(ss,:),'LineWidth',lw)];
            substrs{ss} = sprintf('S%02d',sublist(ss));
        end

        set(gca, 'FontSize', fs, 'XLim',[0 max(tax)])
        ylim(lims);
        plot(get(gca,'XLim'),[0,0],'-','Color',[0.8, 0.8, 0.8]);
        for ee = 1:length(evts2plot)
            plot([evts2plot(ee),evts2plot(ee)],lims,'-','Color',[0.8, 0.8, 0.8]);
        end
        
        if vi==1
            xlabel('Time(s)');
            ylabel('Rand - Pred diff');
        end

        if vi==numel(spat_areas)
            legend(lh,substrs,'FontSize', fs);
        end
       
        title(sprintf('%s',ROI_names{spat_areas(vi)}));
      
    end
    set(gcf,'Color','w')
    set(gcf,'Position',[200,200,1800,1200] );
end
