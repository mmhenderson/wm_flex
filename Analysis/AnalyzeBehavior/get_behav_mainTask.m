% Plot performance on main task, comparing two conditions
% Accuracy and reaction time.
% Also performs stats comparing two conditions.

% Used to create Figure 1B.
%%
close all
clear

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
figpath = fullfile(exp_path,'figs');

sublist = [2,3,4,5,6,7];
nSubj = length(sublist);
% condition colors
col = [125, 93, 175; 15, 127, 98]./255;

condlabs = {'Predictable','Random'};

nCond = 2;
nSess = 2;
nRuns = 10;
nTrials = 20;

% save average stats for each sub and condition
acc_each_run = nan(nSubj, nCond, nSess, nRuns);
RT_each_trial = nan(nSubj, nCond, nSess, nRuns, nTrials);
acc = nan(nSubj, nCond);
RT_mean = nan(nSubj, nCond);

%% load data for all subjects
for si=1:length(sublist)   
    
    ss=sublist(si);
    
    allresp = [];
    allcorrectresp = [];
    alltask = [];
    allrts = [];
   
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
            
        end
    end
    
    % compute condition-specific stats, across all runs this subject.
    
    for cc=1:2
        these_resp = allresp(alltask==cc);
        these_correct_resp =allcorrectresp(alltask==cc);
        these_rts = allrts(alltask==cc);
        
        acc(si,cc) = nanmean(these_resp==these_correct_resp);
        RT_mean(si,cc) = nanmean(these_rts);
        
    end
end

%% plot accuracy overall on two tasks (bar plot)
gray_colors = gray(nSubj+1);
sub_line_color = gray_colors(5,:);

figure;hold all;
meanvals = mean(acc,1);
if nSubj>1
    sevals = std(acc,[],1)./sqrt(nSubj);
else
    sevals = nan(size(meanvals));
end

% first make the actual bar plot
b = bar(1:nCond, meanvals);
b.FaceColor='flat';
b.EdgeColor='flat';
for cc=1:nCond
    b.CData(cc,:) = col(cc,:);
    errorbar(b.XData(cc),meanvals(cc),sevals(cc),'Marker','none',...
            'LineStyle','none','LineWidth',1.5,'Color',col(cc,:),'CapSize',0);
end

% now adding single subjects in gray underneath
for ss = 1:nSubj

    h= plot(1:nCond, acc(ss,:),'-','Color',sub_line_color,'LineWidth',1.5); 
    uistack(h,'bottom');
end
ylim([0.5, 1]);
xlim([0,nCond+1]);
set(gca,'XTick',[1,2],'XTickLabels',condlabs);

title('Accuracy over all runs');
ylabel('Accuracy');
set(gcf,'Color','w');
set(gcf,'Position',[400,400,400,600]);
saveas(gcf,fullfile(figpath,'Acc_allsubs_bars.pdf'),'pdf');
%% compare accuracy across conditions

realvals = acc;
[h,p_ttest,ci,stats] = ttest(realvals(:,1),realvals(:,2));

%% print out accuracy results
fprintf('\n');
for cc=1:nCond
    fprintf('Accuracy for %s: %.2f +/- %.2f\n', condlabs{cc},mean(acc(:,cc))*100, std(acc(:,cc))/sqrt(nSubj)*100)
end
fprintf('parametric paired test tstat=%.3f, df=%d, p value: %.3f\n',stats.tstat,stats.df,p_ttest);
fprintf('num subj showing same effect: %d\n',sum(acc(:,1)>acc(:,2)));

%% plot RT overall on two tasks (bar plot)

gray_colors = gray(nSubj+1);
sub_line_color = gray_colors(5,:);

figure;hold all;
meanvals = mean(RT_mean,1);
if nSubj>1
    sevals = std(RT_mean,[],1)./sqrt(nSubj);
else
    sevals = nan(size(meanvals));
end

% first make the actual bar plot
b = bar(1:nCond, meanvals);
b.FaceColor='flat';
b.EdgeColor='flat';
for cc=1:nCond
    b.CData(cc,:) = col(cc,:);
    errorbar(b.XData(cc),meanvals(cc),sevals(cc),'Marker','none',...
            'LineStyle','none','LineWidth',1.5,'Color',col(cc,:),'CapSize',0);
end

% now adding single subjects in gray underneath
for ss = 1:nSubj

    h= plot(1:nCond, RT_mean(ss,:),'-','Color',sub_line_color,'LineWidth',1.5); 
    uistack(h,'bottom');
end
ylim([0, 1.5]);
xlim([0,nCond+1]);
set(gca,'XTick',[1,2],'XTickLabels',condlabs);

title('RT over all runs');
ylabel('RT (s)');
set(gcf,'Color','w');
set(gcf,'Position',[400,400,400,600]);
saveas(gcf,fullfile(figpath,'RT_allsubs_bars.pdf'),'pdf');

%% compare RT across conditions

realvals = RT_mean;
[h,p_ttest,ci,stats] = ttest(realvals(:,1),realvals(:,2));

%% print out RT results
fprintf('\n');
for cc=1:nCond
    fprintf('RT for %s: %.2f +/- %.2f sec\n', condlabs{cc},mean(RT_mean(:,cc)), std(RT_mean(:,cc))/sqrt(nSubj))
end
fprintf('parametric paired test tstat=%.3f, df=%d, p value: %.3f\n',stats.tstat,stats.df,p_ttest);
fprintf('num subj showing same effect: %d\n',sum(RT_mean(:,1)<RT_mean(:,2)));
