% measure performance on motor cortex button press task
clear
close all

subinit_big = {'BX','BR','CI','CA','CH','AV','CP'};
subnum_big = [1,2,3,4,5,6,7];
sub2do=[2,3,4,5,6,7];
nSubj=length(sub2do);

% find my root directory 
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
root = mypath(1:filesepinds(end-nDirsUp+1));

nRunsMax=8;
allacc = nan(nSubj,nRunsMax);
allrt = nan(nSubj, nRunsMax);
for si = 1:numel(sub2do)
    
    ss=subnum_big(sub2do(si));
    % get subject information
    substr = sprintf('S%02d',subnum_big(ss));
    ii=0;
    behav_folders = dir(fullfile(root,'DataBehavior',substr));
    for se=1:length(behav_folders)
        behav_files = dir(fullfile(behav_folders(se).folder,behav_folders(se).name,'*DigitLoc*.mat'));
        if ~isempty(behav_files)
            load( fullfile(behav_files(1).folder, behav_files(1).name))
            for rr=1:length(TheData)
                resp=TheData(rr).data.Response;
                if ss==6 && rr==1
                    % this is a run where subject reported pressing 2 instead of 4, but
                    % was otherwise good at task
                    resp(resp==2) = 4;
                end
                resp(resp==4) = 2;
                acc = mean(resp==TheData(rr).p.WhichFinger);
                ii=ii+1;
                allacc(si,ii) = acc;
                allrt(si,ii) = nanmean(TheData(rr).data.RespTimeFromOnset);
            end
        end
    end
end
   
avg_each_sub = nanmean(allacc,2);
rt_each_sub = nanmean(allrt,2);
assert(~any(isnan(avg_each_sub)))

fprintf('mean performance %.3f +/- %.3f\n',mean(avg_each_sub),std(avg_each_sub)/sqrt(nSubj))
fprintf('mean RT (s) %.3f +/- %.3f\n',mean(rt_each_sub),std(rt_each_sub)/sqrt(nSubj))