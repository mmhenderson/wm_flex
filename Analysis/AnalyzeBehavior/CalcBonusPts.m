% Calculate subjects' bonus payment compensation - print out the number of
% extra dollars to give per session. 

%% 

close all
clear

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
root = mypath(1:filesepinds(end-nDirsUp+1));

sublist = [2:7];
nSubj = length(sublist);
condlabs = {'Predictable','Random'};

nCond = 2;
nSess = 2;
nRuns = 10;

% how many dollars per point? 
dollars_per_pt = 1.00;

% how many bonus trials occured in every session and condition?
nBonusPerSessPerCond = 10;

% save average stats for each sub and condition
totalBonusCorrect = nan(nSubj, nSess, nCond);
totalBonusIncorrect = nan(nSubj, nSess, nCond);

calc_dollars_per_sess = nan(nSubj, nSess);

%% load data for all subjects
for si=1:numel(sublist)
    ss=sublist(si);
    
    for se=1:nSess
        
        bonus_correct = zeros(nCond, 1);
        bonus_total = zeros(nCond, 1);
    
        if ss==6 && se==2
            files = dir(sprintf('%sDataBehavior/S0%d/Session%d/*MainTask*.mat',root,ss,3));
        else
            files = dir(sprintf('%sDataBehavior/S0%d/Session%d/*MainTask*.mat',root,ss,se));
        end
        load([files(1).folder '/' files(1).name])
        
        for rr=1:nRuns
            
            p = TheData(rr).p;
            t = TheData(rr).t;
            data = TheData(rr).data;
            
            for cc = 1:nCond
                
                inds = p.Task==cc & p.BonusTrial==1;
                
                nBonusTotal = sum(inds);
                nBonusCorrect = sum(p.CorrectResponse(inds)==data.Response(inds));
                
                bonus_correct(cc,1) = bonus_correct(cc,1) + nBonusCorrect;
                bonus_total(cc,1) = bonus_total(cc,1) + nBonusTotal;
            end                                            
        end

        assert(all(bonus_total==nBonusPerSessPerCond))
        totalBonusCorrect(si,se,:) = bonus_correct; 
        calc_dollars_per_sess(si,se) = sum(totalBonusCorrect(si,se,:))*dollars_per_pt;

    end
    
    %% print result
    fprintf('\nSubject %d:\n',ss)
    for se = 1:nSess
        fprintf('  Session %d: %d bonus pts total (%d %s, %d %s)\n',se, ...
            sum(totalBonusCorrect(si,se,:)), totalBonusCorrect(si,se,1), condlabs{1},totalBonusCorrect(si,se,2), condlabs{2});
        fprintf('     $ %.2f\n', calc_dollars_per_sess(si,se));
    end
end
npersubj=sum(calc_dollars_per_sess,2);
fprintf('average bonus per subj = %.2f +/- %.2f\n',mean(npersubj),std(npersubj))