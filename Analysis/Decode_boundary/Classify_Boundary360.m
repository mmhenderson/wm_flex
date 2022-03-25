%% Classifying the orientation of the disk boundary. 
% Treating the orientations as a full 360 deg space, so 0/180 are not same
% (because the disk sides are different)
% Average these two decoding results.
% This script saves a mat file with decoding performance, can plot it using
% plotClassResults_Boundary.m
% note this analysis is not in our paper.
%%
clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
curr_dir = pwd;
filesepinds = find(curr_dir==filesep);
nDirsUp = 2;
exp_path = curr_dir(1:filesepinds(end-nDirsUp+1));

nVox2Use = 10000;    % this is the max number of vox to use, so if it's very big we're using all the voxels.

condLabStrs = {'Predictable','Random'};
nConds = length(condLabStrs);

class_str = 'normEucDist';

tr_avg_range = [6, 12];

%% loop over subjects
for ss=1:length(sublist)

    substr = sprintf('S%02d',sublist(ss));
    
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(curr_dir,'Decoding_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('ClassifyBoundary360_%s_%dvox_%s.mat',class_str,nVox2Use,substr));

   
    v2do = 1:length(ROI_names);
    for vv = v2do
%     for vv = 1:length(ROI_names)
        
        %% pull out the data for main task

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<1
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        
        for cc = 1:nConds
            
            % always using label of the boundary that was physically shown
            % during preview (so for random, this is the random preview,
            % not the final resp)
            if cc==1
                boundLabs = mainSig(vv).boundPos;
            else
                boundLabs = mainSig(vv).randBoundPos;
            end
            condLabs = mainSig(vv).condLabs;
            runLabs = mainSig(vv).runLabs;
            % take out just condition of interest
            trials2use = condLabs==cc;
            
            boundLabs = boundLabs(trials2use);
            boundLabs = round(boundLabs,1);
%             % only want boundary ORIENTATION here, so 90/270 are same
%             boundLabs = round(mod(boundLabs,180),1);
            mainDat = mainSig(vv).dat_by_TR;
            % make sure to add one to the tr index range, because the data
            % here is zero-indexed (the first TR is the TR at which event
            % started, so time zero)
            mainDat = mean(mainDat(:,tr_avg_range(1)+1:tr_avg_range(2)+1,:),2);
            mainDat = mainDat(trials2use,:);

            % bin these for classifier - want 4 bins that are roughly centered at
            % 0, 45, 90, 135.
            binLabs = zeros(size(boundLabs));
            nbins=8;
            bin_centers=0:45:360;

            bin_size=diff(bin_centers(1:2));
            for bb=1:nbins
                inds_this_bin = abs(boundLabs-(bin_centers(bb)-0.0001))<bin_size/2 | ...
                                abs((boundLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
                binLabs(inds_this_bin) = bb;
            end
            assert(~any(binLabs==0))
            
            % subtract mean over voxels 
            mainDat = mainDat - repmat(mean(mainDat,2), 1, size(mainDat, 2));

            if vv==v2do(1) && cc==1
                % preallocate array here
                % ngroups is the two binary classifications that get done -
                % 0 versus 90, and 45 versus 135
                nGroups=4;
                allacc = nan(length(ROI_names), nConds, nGroups);
                alld = nan(length(ROI_names), nConds, nGroups);
            end
            
            % doing cross-validation over sessions here, becuase this keeps
            % the training sets perfectly balanced.
            sessLabs = ones(size(runLabs));
            sessLabs(runLabs>10) = 2;
            
            cvLabs = sessLabs(trials2use);
            nCV = numel(unique(cvLabs));

            dat2use = mainDat;

            %% voxel selection from each training set 
         
            if ~isempty(nVox2Use) && nVox2Use<size(dat2use,2)
                numcores = 8;
                if isempty(gcp('nocreate'))
                    parpool(numcores);
                end
                fprintf('running voxel selection f-test for %s %s - %s condition\n',substr, ROI_names{vv}, condLabStrs{cc})
                voxStatTable = zeros(size(dat2use,2),nCV);
                for rr = 1:nCV
                    inds = cvLabs~=rr;
                    pvals = zeros(size(dat2use,2), 1);
                    dat = dat2use(inds,:);
                    lab = binLabs(inds,:);
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

            % going to do four different classifications - bin 1 versus bin 3,
            % and bin 2 versus bin 4, etc
            groups = [1,3; 2,4; 5,7; 6,8];
            for xx = 1:nGroups

                inds2use = binLabs==groups(xx,1) | binLabs==groups(xx,2);

                %% define train and test set 

                trnDat = dat2use(inds2use,:);
                trnLabs = binLabs(inds2use,:);
                trnCV = cvLabs(inds2use,:);

                tstDat = dat2use(inds2use,:);
                tstLabs = binLabs(inds2use,:);
                tstCV = cvLabs(inds2use,:);

                % run classifier using custom code to cross-validate
                [~,~,predLabs] = my_classifier_cross(trnDat,trnLabs,...
                    trnCV,tstDat, tstLabs,...
                    tstCV,class_str,100,nVox2Use_now,voxStatTable,1);

                acc = mean(predLabs==tstLabs);
                dprime = get_dprime(predLabs, tstLabs,tstLabs);

                allacc(vv,cc,xx) = acc;
                alld(vv,cc,xx) = dprime;
                
            end
        end
        

    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allacc','alld');

end