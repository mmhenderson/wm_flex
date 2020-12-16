%% plot reconstructions
% training on spatial working memory (SWM) mapping task
% testing on main task conditions
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

% names of the ROIs - note the motor areas have been changed to nicer names
% here for plotting.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'S1','M1','PMc',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA','M1/S1 all'};

plot_order1 = [1:5,10,11,6:9,12:14]; 
vismotor_names = ROI_names(plot_order1);
plot_order_all = [plot_order1];
vismotor_inds = find(ismember(plot_order_all,plot_order1));

nROIs = length(plot_order_all);

plotVisMotorFids = 1;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

nVOIs = length(plot_order_all);
trDur = 0.80;

col = [125, 93, 175; 15, 127, 98]./255;

ylims_fid = [-0.1, 0.25];
ylims_bias=[-90,90];
lw = 1;
fs=20;



evts2plot = [3.5, 4.5, 16.5, 18.5];

sig_heights = [0.19, 0.21, 0.23];
diff_col=[0.5, 0.5, 0.5];

alpha_vals=[0.05,0.01,0.001];
alpha=alpha_vals(1);
alpha_ms = [8,16,24];

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
        
    end
  
    condlabs = allchanresp(1).condLabs;
   
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
               
           end
       end
    end  
end

% compute some basic stats here
vals = fidelity;
meanvals = squeeze(mean(vals,1));
semvals = squeeze(std(vals,[],1)./sqrt(nSubj));
%% pairwise condition comparisons
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 143545;
rng(rndseed,'twister')
nPermIter=1000;
real_sr_stat = nan(nROIs,nTRs_out);
rand_sr_stat = nan(nROIs, nTRs_out, nPermIter);

for vv=1:nROIs
    for tr=1:nTRs_out
        realvals = squeeze(vals(:,vv,:,tr));

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
