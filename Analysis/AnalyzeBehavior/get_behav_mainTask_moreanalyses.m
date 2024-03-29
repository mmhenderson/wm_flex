% Plot performance on main task, both conditions

% Includes additional behavioral analyses not reported in our paper.
%%
close all
clear

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

addpath(fullfile(exp_path,'Analysis','plotting_utils'))
figpath = fullfile(exp_path,'figs');

sublist = [2,3,4,5,6,7];
nSubj = length(sublist);
subcolors = viridis(nSubj+1);
% condition colors
col = [125, 93, 175; 15, 127, 98]./255;


condlabs = {'Predictable','Random'};

nCond = 2;
nSess = 2;
nRuns = 10;
nTrials = 20;
ylimhist_rt = [0,100];
nhistbins = 20;

% save average stats for each sub and condition
acc_each_run = nan(nSubj, nCond, nSess, nRuns);
RT_each_trial = nan(nSubj, nCond, nSess, nRuns, nTrials);
acc = nan(nSubj, nCond);
RT_mean = nan(nSubj, nCond);
RT_stdev = nan(nSubj, nCond);

% look at performance in different "difficulty" bins, meaning the distance
% between the target and the boundary.
dist_bin_size = 10;
nBins = 90/dist_bin_size;   
bin_edges = dist_bin_size:dist_bin_size:90;
%  note this space is only 90 deg wide, because the target can never be 
% more than 90 degrees from the bound (abs)

acc_binned = nan(nSubj, nCond, nBins);
RT_mean_binned = nan(nSubj, nCond, nBins);
RT_stdev_binned = nan(nSubj, nCond, nBins);

nperbin = zeros(nSubj, nCond, nBins);
%% load data for all subjects
for si=1:length(sublist)   
    
    ss=sublist(si);
    
    allresp = [];
    allcorrectresp = [];
    alltask = [];
    allrts = [];
    alltarg = [];
    allbound = [];
    
    for se=1:nSess
        

        files = dir(sprintf('%sDataBehavior/S0%d/Session%d/*MainTask*.mat',exp_path,ss,se));
        if ss==6 && se==2
            files = dir(sprintf('%sDataBehavior/S0%d/Session3/*MainTask*.mat',exp_path,ss));
        end
        if isempty(files)
            fprintf('no data for sub %d sess %d\n',ss,se);
            continue
        end
        load([files(1).folder '/' files(1).name])
        
        for rr=1:nRuns
            
            if rr>length(TheData)
                fprintf('no data for sub %d sess %d run %d\n',ss,se,rr);
                continue
            end
            p = TheData(rr).p;
            t = TheData(rr).t;
            data = TheData(rr).data;
            
            for cc = 1:nCond
                
                % get stats about this run, this condition
                these_resp = data.Response(p.Task==cc);
                these_correct_resp = p.CorrectResponse(p.Task==cc);
                acc_each_run(si,cc,se,rr) = nanmean(these_resp==these_correct_resp); 
                these_rts = t.respTimeFromStart(p.Task==cc);
                RT_each_trial(si,cc,se,rr,1:length(these_rts)) = t.respTimeFromStart(p.Task==cc);
            end
            
            % make a long list of stuff to concat across all runs
            allresp = [allresp; data.Response];
            allcorrectresp = [allcorrectresp; p.CorrectResponse];
            alltask = [alltask; p.Task];
            allrts = [allrts; t.respTimeFromStart];
            alltarg = [alltarg; p.TargPosition];
            allbound = [allbound; p.BoundPosition];
            
        end
    end
    
    % compute condition-specific stats, across all runs this subject.
    
    for cc=1:2
        these_resp = allresp(alltask==cc);
        these_correct_resp =allcorrectresp(alltask==cc);
        these_rts = allrts(alltask==cc);
        
        acc(si,cc) = nanmean(these_resp==these_correct_resp);
        RT_mean(si,cc) = nanmean(these_rts);
        RT_stdev(si,cc) = nanstd(these_rts);
        
        % define absolute distance between target and boundary.
        dist1 = abs(alltarg(alltask==cc)-allbound(alltask==cc));
        dist2 = abs(180-dist1);
        dist3 = abs(360-dist1);
        mindist = min([dist1, dist2, dist3],[],2);
        % now looking within difficulty/distance bins
        for bb = 1:nBins
          
           if bb>1
               this_bin_inds = mindist<=bin_edges(bb) & mindist>bin_edges(bb-1);
           else
               this_bin_inds = mindist<=bin_edges(bb);
           end
           nperbin(ss,cc,bb) = sum(this_bin_inds);
           
           binresp = these_resp(this_bin_inds);
           bincorrectresp = these_correct_resp(this_bin_inds);
           acc_binned(si,cc,bb) = nanmean(binresp==bincorrectresp);
           RT_mean_binned(si,cc,bb) = nanmean(these_rts(this_bin_inds));
           RT_stdev_binned(si,cc,bb) = nanstd(these_rts(this_bin_inds));
           
        end
    end
end


%% plot RT histograms, each condition
figure;hold all;

for cc = 1:nCond
    
    subplot(2,1,cc);hold all;
    
    dat = squeeze(RT_each_trial(:,cc,:,:,:));
    dat = dat(~isnan(dat));
    % print the minimum values
    disp(min(dat))
    histogram(dat,200,'FaceColor',col(cc,:),'EdgeColor',col(cc,:));
    title(condlabs{cc})
    plot([0.2, 0.2],get(gca,'YLim'),'k')
    xlabel('RT (sec)')
    ylabel('Number of trials')
    xticks(0:0.2:3)
    xlim([0,3])
    set(gca,'TickDir','out')
end
saveas(gcf,fullfile(figpath,'RT_hist.pdf'),'pdf');
    
%% plot accuracy as a function of run number
figure;hold all;

for cc = 1:nCond
    ll = cell(nSubj,1);
    lh = [];
    subplot(2,1,cc);hold all;
    for ss = 1:nSubj    
        ll{ss} = sprintf('S%02d',ss);
        dat = acc_each_run(ss,cc, :,:); % this is [nSess x nRuns];
        dat = squeeze(dat)';    % this is [nRuns x nSess];
        dat = dat(:);   % this is [nRunsTotal x 1]
        h= plot(1:nRuns*nSess, dat,'-','Color',subcolors(ss,:)); 
        lh = [lh, h];
        
    end
    title(condlabs{cc});
    ylim([0,1.1]);
    xlim([0,nRuns*nSess+1]);
    xlabel('Runs in whole expt')
    plot(get(gca,'XLim'), [0.5, 0.5],'--','Color',[0.8 0.8, 0.8]);
    legend(lh,ll);
end
suptitle('Accuracy over time');

set(gcf,'Color','w');

%% plot accuracy based on difficulty (distance from stim to boundary)
figure; hold all;
for cc = 1:nCond
    subplot(2,1,cc);hold all;
    ll = cell(nSubj,1);
    lh = [];
    for ss = 1:nSubj
        ll{ss} = sprintf('S%02d',ss);
        dat = squeeze(acc_binned(ss,cc,:)); %[nBins] values

        h= plot(1:nBins, dat,'-','Color',subcolors(ss,:)); 
        lh = [lh, h];
        
    end
    title(condlabs{cc});
    ylim([0,1.1]);
    xlim([0,nBins+1]);
    xlabel('Angular difference')
    set(gca,'XTick',1:nBins,'XTickLabels',bin_edges-dist_bin_size/2);
    plot(get(gca,'XLim'), [0.5, 0.5],'--','Color',[0.8 0.8, 0.8]);
    legend(lh,ll);
end

suptitle('Accuracy versus difficulty');
set(gcf,'Color','w');

%% plot RT as a function of run number
figure;hold all;

for cc = 1:nCond
    ll = cell(nSubj,1);
    lh = [];
    subplot(2,1,cc);hold all;
    ylim([0,3.1]);
    for ss = 1:nSubj   
        ll{ss} = sprintf('S%02d',ss);
        curr_x = 0;
        for se  =1:nSess
            for rr = 1:nRuns
        
                dat = squeeze(RT_each_trial(ss,cc, se,rr,:));
                dat = dat(~isnan(dat));
                
                h= plot(curr_x+1:curr_x+length(dat), dat,'o','Color',subcolors(ss,:)); 
                if se==1 && rr==1
                    lh = [lh, h];
                end
                curr_x = curr_x+length(dat);
%                 line([curr_x+0.5,curr_x+0.5], get(gca,'YLim'),'Color',[0.8, 0.8, 0.8]); 
            end
        end
        
    end
    title(condlabs{cc});    
    xlim([0,curr_x+1]);
    xlabel('Trials in whole expt');
    legend(lh,ll);
end
suptitle('RT over time');
set(gcf,'Color','w');


%% plot RT based on difficulty (distance from stim to boundary)
figure; hold all;
for cc = 1:nCond
    subplot(2,1,cc);hold all;
    ll = cell(nSubj,1);
    lh = [];
    for ss = 1:nSubj
        ll{ss} = sprintf('S%02d',ss);
        dat = squeeze(RT_mean_binned(ss,cc,:)); %[nBins] values
        errdat = squeeze(RT_stdev_binned(ss,cc,:));
        h = errorbar(1:nBins, dat, errdat, 'Color',subcolors(ss,:));
%         h= plot(1:nBins, dat,'-','Color',subcolors(ss,:)); 
        lh = [lh, h];
        
    end
    title(condlabs{cc});
    ylim([0,3.2]);
    xlim([0,nBins+1]);
    xlabel('Angular difference')
    ylabel('RT (sec)');
    set(gca,'XTick',1:nBins,'XTickLabels',bin_edges-dist_bin_size/2);
    legend(lh,ll);
end

suptitle('RT versus difficulty');
set(gcf,'Color','w');