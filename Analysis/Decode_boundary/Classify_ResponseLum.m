% MMH 3/12/20
% classifying which luminance was the side of the response (note this is
% opposite finger for different mappings)

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
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyResponseLum_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

   
    v2do = 1:length(ROI_names);
    for vv = v2do
%     for vv = 1:length(ROI_names)
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        for cc = 1:nConds
            
%             respLabs = mainSig(vv).RespActual;
            respLabs = mainSig(vv).CorrectResp;
            
            condLabs = mainSig(vv).condLabs;
            runLabs = mainSig(vv).runLabs;
            % getting rid of any trials with no response here
            % also taking out just the relevant condition!!
            trials2use = ~isnan(respLabs) & respLabs~=0 & condLabs==cc;
            
            respLabs = respLabs(trials2use);
            runLabs = runLabs(trials2use);
            sessLabs = ones(size(runLabs));
            sessLabs(runLabs>10) = 2;
            cvLabs = sessLabs;
            
            % classifying response luminance - which luminance corresponded
            % to the correct answer?
            lumLabs = respLabs;
            lumLabs(sessLabs==2) = 3-lumLabs(sessLabs==2);
            
            mainDat = mainSig(vv).dat_avg;
            mainDat = mainDat(trials2use,:);

            % subtract mean over voxels 
%             mainDat = mainDat - repmat(mean(mainDat,2), 1, size(mainDat, 2));

            if vv==v2do(1) && cc==1
                % preallocate array here
                allacc = nan(length(ROI_names), nConds);
                alld = nan(length(ROI_names), nConds);
            end

            nCV = numel(unique(cvLabs));

            dat2use = mainDat;

            %% voxel selection from each training set 
            % for this voxel selection I'm using trials from all conditions, but
            % leaving out one session at a time. This gives a list of voxels to use
            % for each fold of cross validation. The same voxels are used
            % regardless of which condition we are using for classification. Think
            % this makes the condition comparisons more fair. Also saves time
            % because we only need to run this once.
            if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                voxStatTable = zeros(size(dat2use,2),nCV);
                for rr = 1:nCV
                    inds = cvLabs~=rr;
                    pvals = zeros(size(dat2use,2), 1);
                    dat = dat2use(inds,:);
                    lab = lumLabs(inds,:);
                    parfor vx = 1:size(dat2use,2)
                         % choose the voxels        
                       [pvalue, stats] = anovan(dat(:,vx), lab,'display','off');
                       pvals(vx) = pvalue;
                    end 
                    voxStatTable(:,rr) = pvals;
                end
                nVox2Use_now = nVox2Use;
            else            
                % put in a placeholder here because using all voxels
                voxStatTable = zeros(size(dat2use,2),nCV);
                nVox2Use_now = [];
            end

            %% define train and test set 

            % same data here because we're not cross-generalizing or anything
            trnDat = dat2use;
            trnLabs = lumLabs;
            trnCV = cvLabs;

            tstDat = dat2use;
            tstLabs = lumLabs;
            tstCV = cvLabs;

            %% run the classifier w/ balancing if needed


            [~,~,predLabs] = my_classifier_cross(trnDat,trnLabs,...
                trnCV,tstDat, tstLabs,...
                tstCV,class_str,100,nVox2Use_now,voxStatTable,1);

            acc = mean(predLabs==tstLabs);
            dprime = get_dprime(predLabs, tstLabs,tstLabs);

            allacc(vv,cc) = acc;
            alld(vv,cc) = dprime;
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld');

end