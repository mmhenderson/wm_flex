% measure performance on wedge mapping task
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

nRunsMax=16;
allhr = nan(nSubj,nRunsMax);

for si = 1:numel(sub2do)
    
    ss=subnum_big(sub2do(si));
    % get subject information
    substr = sprintf('S%02d',subnum_big(ss));
    ii=0;
    behav_folders = dir(fullfile(root,'DataBehavior',substr));
    for se=1:length(behav_folders)
        behav_files = dir(fullfile(behav_folders(se).folder,behav_folders(se).name,'*sIEM_1D*.mat'));
        for rr=1:numel(behav_files)
            load( fullfile(behav_files(rr).folder, behav_files(rr).name))
                     
            ii=ii+1;
            allhr(si,ii) = mean(p.runHits);
            
        end
    end
    clear TheData
end
   
avg_each_sub = nanmean(allhr,2);
assert(~any(isnan(avg_each_sub)))

fprintf('mean hit rate %.3f +/- %.3f\n',mean(avg_each_sub),std(avg_each_sub)/sqrt(nSubj))
