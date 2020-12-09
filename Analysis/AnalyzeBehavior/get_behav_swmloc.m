% measure performance on spatial working memory localizer task
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

nRunsMax=10;
allerr = nan(nSubj,nRunsMax);

for si = 1:numel(sub2do)
    
    ss=subnum_big(sub2do(si));
    % get subject information
    substr = sprintf('S%02d',subnum_big(ss));
    ii=0;
    behav_folders = dir(fullfile(root,'DataBehavior',substr));
    for se=1:length(behav_folders)
        behav_files = dir(fullfile(behav_folders(se).folder,behav_folders(se).name,'*SWMLoc*.mat'));
        if ~isempty(behav_files) && ~contains(behav_files(1).name,'TRAINING')
%         if ~isempty(behav_files) 
            load( fullfile(behav_files(1).folder, behav_files(1).name))
            for rr=1:length(TheData)               
                ii=ii+1;
                allerr(si,ii) = mean(TheData(rr).data.Error);
            end
        end
    end
    clear TheData
end
   
avg_each_sub = nanmean(allerr,2);
assert(~any(isnan(avg_each_sub)))

fprintf('mean error %.1f +/- %.1f deg\n',mean(avg_each_sub),std(avg_each_sub)/sqrt(nSubj))
