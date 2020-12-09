% measure performance on retinotopic mapping (rotating wedge) task for my
% subjects in this experiment
clear
close all

subinit_big = {'BX','BR','CI','CA','CH','AV','CP'};
subnum_big = [1,2,3,4,5,6,7];
sub2do=[2,3,4,5,6,7];
nSubj=length(sub2do);

% find my root directory 
ret_path = '/mnt/neurocube/local/serenceslab/Doreti/';

nRunsMax=10;
allacc = nan(nSubj,nRunsMax);
for si = 1:numel(sub2do)
    
    ss=subnum_big(sub2do(si));
    % get subject information
    subinit = subinit_big{ss};
    substr = sprintf('S%02d',subnum_big(ss));
    
    behav_folder = fullfile(ret_path,'RAW',subinit,'behav');
    behav_files = dir(fullfile(behav_folder,'*wedge*.mat'));
    for rr=1:numel(behav_files)
       fn2load = fullfile(behav_files(rr).folder, behav_files(rr).name);
       load(fn2load);
       allacc(si,rr) = p.acc;
    end
end
   
avg_each_sub = nanmean(allacc,2);
assert(~any(isnan(avg_each_sub)))

fprintf('mean performance %.2f +/- %.2f\n',mean(avg_each_sub),std(avg_each_sub)/sqrt(nSubj))