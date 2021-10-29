% MMH 3/12/20
% classifying RESPONSE - WHICH FINGER DID THEY EVENTUALLY PRESS??
    % note this is counter-balanced w/r/t COLOR of the side.

clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;
filesepinds = find(curr_dir==filesep);
nDirsUp = 2;
exp_path = curr_dir(1:filesepinds(end-nDirsUp+1));

nVox2Use = 10000;

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

class_str = 'normEucDist';
% class_str = 'svmtrain_lin';

dbstop if error
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    fn2load = fullfile(exp_path,'Samples',sprintf('DWMLocSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponse_TrnDWMLoc_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

   
    v2do = 1:length(ROI_names);
    for vv = v2do
%     for vv = 1:length(ROI_names)
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        for cc = 1:nConds
            
            %% first organize my training data (digit localizer)
            
            % using all trials because accuracy is generally very high. Can
            % also remove the incorrect trials, but that would un-balance
            % the training set.
            trials2use = ones(size(locSig(1).ActDigit,1),1)==1;
%             trials2use = locSig(vv).ActDigit == locSig(vv).ExpDigit;
            trnLabs = locSig(vv).ExpDigit(trials2use);
            trnRuns = locSig(vv).runLabs(trials2use);
            trnDat = locSig(vv).dat_avg_targ(trials2use,:);
%             trnDat = trnDat - repmat(mean(trnDat,2), 1, size(trnDat,2));
            
            %% now organize my testing data (main task)
            
            tstLabs = mainSig(vv).RespActual;
            condLabsTest = mainSig(vv).condLabs;
            tstRuns = mainSig(vv).runLabs;
            % getting rid of any trials with no response here
            % also taking out just the relevant condition!!
            trials2use = ~isnan(tstLabs) & tstLabs~=0 & condLabsTest==cc;
            
            tstLabs = tstLabs(trials2use);

            mainDat = mainSig(vv).dat_avg;
            mainDat = mainDat(trials2use,:);

            % subtract mean over voxels 
%             mainDat = mainDat - repmat(mean(mainDat,2), 1, size(mainDat, 2));

            tstRuns = tstRuns(trials2use);
            nRuns = numel(unique(tstRuns));

            if vv==v2do(1) && cc==1
                % preallocate array here
                allacc = zeros(length(ROI_names), nConds);
                alld = zeros(length(ROI_names), nConds);
            end
            
            
            tstDat = mainDat;

            %% voxel selection from each training set 
            % usng the training set only, find voxels that are selective
            % for response (which finger?)
            if ~isempty(nVox2Use) && nVox2Use<size(trnDat,2)
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                
                pvals = zeros(size(trnDat,2), 1);
                dat = trnDat;
                lab = trnLabs;
                parfor vx = 1:size(trnDat,2)
                     % choose the voxels        
                   [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                   pvals(vx) = pvalue;
                end
                assert(~any(isnan(pvals)))
                % use the correct column of pTable to sort
                [~, ind] = sort(pvals, 'ascend'); 
                voxelindsuse = ind(1:nVox2Use);
            else            
                % put in a placeholder here because using all voxels
                voxelindsuse = 1:size(trnDat,2);
            end

            %% run the classifier w/ balancing if needed

            [predLabs,~] = normEucDistClass(trnDat(:,voxelindsuse),tstDat(:,voxelindsuse),trnLabs);
            
            acc = mean(predLabs==tstLabs);
            dprime = get_dprime(predLabs, tstLabs,tstLabs);

            allacc(vv,cc) = acc;
            alld(vv,cc) = dprime;
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld');

end