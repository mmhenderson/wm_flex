function [acc,dprime,predLabs,distEachClass] = my_classifier_cross_wconf(trndat,trnlabs,trnruns,tstdat,tstlabs,tstruns,classstr,nBalanceIter,nVox2Keep,voxStatTable,resamp)
% classify voxel activation patterns into n groups using a linear method.

% INPUTS: 
% trndat: [nTrials x nVoxels], data for each train trial
% trnlabs: [nTrials x 1], labels for each train trial, can be any unique values.
% trnruns: [nTrials x 1], stores the number of run for each train trial (this is used to
%   crossvalidate, could also be session numbers etc).

% tstdat: [nTrials x nVoxels], data for each test trial
% tstlabs: [nTrials x 1], labels for each test trial, can be any unique values.
% tstruns: [nTrials x 1], stores the number of run for each test trial (this is used to
%   crossvalidate, could also be session numbers etc).

% classstr: a string specifying which type of classifier to use. Defaults
%   to a linear discriminant implemented in classify.m. Look down in the
%   code to see what the options are.
% resamp: integer specifying how to resample, if the trn set is unbalanced.
%   if resamp==1, will downsample larger set
%   if resamp==2, will oversamp (leads to bias so probably don't do this!)
%   if resamp==0, no resampling is done.
% nBalanceIter: (integer), is number of iterations of resampling. Defaults to 100.
% voxStatTable: [nVox x nRuns] used to select a subset of voxels to use on 
%   each cross-validation fold. The columns indicate which cross-validation 
%   fold the statistics correspond to. These values could be:
%   - p-value from an anova (make sure from training set dat only)
%   - 1/t-score of voxel in independent visual localizer.
%   We will always choose voxels by their MINIMUM value in this table.
% nVox2Keep: integer for the number of voxels to use, based on values in
%   voxStatTable.

% OUTPUTS:
% acc: proportion of accurate predictions over whole data set.
% dprime: d' computed over whole data set, see function get_dprime.m
% predLabs: [nTrials x 1], predictions from original data labels
%   predLabs is always double, regardless of input type.
% distEachClass [nTrials x nClasses] confidence of the classifier on single trials, based
% on estimate of distance between test trials and mean of trainign set
% groups

% MMH 7/12/18: this version lets us put in a separate training and testing
% set, as long as all the run numbers exist and correspond in both sets. 
% MMH 11/4/20: this version returns the confidence. Note this only works
% for eucDistClass.m and normEucDistClass.m, don't have a reliable way of
% estimating confidence for other kinds of classifiers yet!
%% set up, check inputs 
if nargin<6
    error('you must supply trndat, tstdat, trnlabs, tstlabs, trnruns, and tstruns')
end
if nargin<11; resamp = []; end
if nargin<10; voxStatTable = []; end
if nargin<9; nVox2Keep = []; end
if nargin<8; nBalanceIter = []; end
if nargin<7; classstr = []; end

unruns=unique(trnruns);
if any(unruns~=unique(tstruns))
    error('you must have all the same unique run numbers in trn and testing sets')
end
nruns=length(unruns);
ntrials_trn = size(trndat,1);
ntrials_tst = size(tstdat,1);
nvox = size(trndat,2);
assert(nvox==size(tstdat,2))

if size(trnlabs,1)~=ntrials_trn || size(trnruns,1)~=ntrials_trn; error('trndat, trnlabs, and trnruns must have same first dimension'); end
if size(tstlabs,1)~=ntrials_tst || size(tstruns,1)~=ntrials_tst; error('tstdat, tstlabs, and tstruns must have same first dimension'); end

if ~isa(trnlabs, 'double')
    trnlabs = double(trnlabs);
end
if ~isa(tstlabs, 'double')
    tstlabs = double(tstlabs);
end

% set default vals
if isempty(nVox2Keep); nVox2Keep = nvox; end
if isempty(resamp);resamp = 1; end
if isempty(classstr); classstr = 'classify_diaglinear';end
if isempty(voxStatTable); voxStatTable = zeros(nvox, nruns); end

if resamp==0
    nBalanceIter = 1;
elseif isempty(nBalanceIter)
    nBalanceIter = 100;
end

if mod(nBalanceIter,1) || nBalanceIter<1; error('nBalanceIter must be an integer value > 0'); end
if nBalanceIter==1; resamp = 0; end
if mod(nVox2Keep, 1); error('nVox2Keep must be an integer value'); end
if nVox2Keep<1 || nVox2Keep>size(trndat,2);  error('nVox2Keep must be in range 1-nVox, or empty'); end
if any(size(voxStatTable)~=[nvox,nruns]); error('voxStatTable must be [nVox x nRuns]'); end
if sum(voxStatTable(:)==0) && nVox2Keep<nvox 
    nVox2Keep = size(trndat,2);
    warning('you specified nVox2Keep as less than the total num voxels, but didn''t specify voxStatTable. Using all voxels.'); 
end

% check to see which svmtrain is on path...by default this would be the
% matlab built in, but better to use libSVM.
if contains(classstr,'svm')
    svm_version = which('svmtrain');
    if ~strcmp(svm_version(end-14:end), 'svmtrain.mexa64')
        error('Verify that you are using the libSVM (3.1) version of "svmtrain". If you think you are using the correct version, you can comment this out and proceed.')
    end
end
if strcmp(classstr, 'eucDist') && exist('eucDistClass','file')~=2
    error('make sure the file eucDistClass.m is on your path.\n%s','A copy of this file is in maggie/mFiles/Classifiers/')
end
if strcmp(classstr, 'normEucDist') && exist('normEucDistClass','file')~=2
    error('make sure the file normEucDistClass.m is on your path.\n%s','A copy of this file is in maggie/mFiles/Classifiers/')
end
if strcmp(classstr, 'ideal_observer') && exist('ideal_observer','file')~=2
    error('make sure the file ideal_observer.m is on your path.\n%s','A copy of this file is in maggie/mFiles/Classifiers/')
end
     
% these are all the groups in the full data set.
balgroups = sort(unique(trnlabs));

% array to store all predicted labels
allpredlabs = zeros(length(tstruns), nBalanceIter);
alldist =zeros(length(tstruns), length(balgroups), nBalanceIter);

%% loop over runs/folds
for cv=1:nruns    
    
    % get my voxels of interest for this cross-validation fold
    if nVox2Keep<nvox
        % use the correct column of pTable to sort
        [~, ind] = sort(voxStatTable(:,cv), 'ascend'); 
        voxelindsuse = ind(1:nVox2Keep);
    else
        voxelindsuse = 1:nvox;
    end

    % get my training and testing sets
    % remember that we're indexing into different arrays here to get the
    % train and test sets. But we still want to cross-validate across runs,
    % so we're only taking a portion of each array at a time.
    tstinds=tstruns==unruns(cv);
    trninds=trnruns~=unruns(cv);
    
    trndat_thisloop = trndat(trninds,voxelindsuse);
    tstdat_thisloop = tstdat(tstinds,voxelindsuse);
    trnlabs_thisloop = trnlabs(trninds);
    tstlabs_thisloop = tstlabs(tstinds);
   
    %see if there's an unbalanced set     
    un=unique(trnlabs_thisloop);
    numeach = sum(repmat(trnlabs_thisloop,1,numel(un))==repmat(un',length(trnlabs_thisloop),1));
   
    %check if any groups are missing entirely
    if numel(balgroups)~=numel(un) || any(balgroups~=sort(un)) 
        error('you are missing trials from one or more label groups, on at least one cross-validation fold. Check your labels!');
    end

    nMax=max(numeach);
    nMin=min(numeach);

    % decide whether we need to balance trn set now
    doBalance = 0;
    if any(numeach~=nMax)
        if resamp>0
            doBalance=1;
        else
            fprintf('Fold %d: training set is unbalanced, but you have opted not to resample. proceeding...\n',cv)
        end
    else
%         fprintf('Fold %d: training set is balanced. Good job!\n', cv)
    end
    
    
    if doBalance
        
        % loop over balancing iterations here.
        labs_tmp = zeros(length(tstlabs_thisloop),nBalanceIter);
        dist_tmp = zeros(length(tstlabs_thisloop),numel(un), nBalanceIter);
        fprintf('Fold %d: training set is unbalanced, starting resampling over %d iters...\n',cv,nBalanceIter);
        parfor bb=1:nBalanceIter
            
            useinds=[];
            if resamp==1
                % want to downsample. Get one possible balanced training set               
                for ii=1:length(un)
                    theseinds=find(trnlabs_thisloop==un(ii));
                    if numel(theseinds)>nMin
                        sampinds = datasample(theseinds,nMin,'Replace',false);
                        useinds = [useinds;sampinds];
                    else
                        useinds = [useinds;theseinds];                    
                    end
                end
                if length(useinds)~=nMin*length(un)
                    error('mistake in resampling')
                end
            elseif resamp==2
                % want to oversample. Get one possible balanced training set                
                for ii=1:length(un)
                    %inds of the actual data
                    theseinds=find(trnlabs_thisloop==un(ii));
                    nOverSamp = nMax-numeach(ii);
                    %inds of the extra data, if any, to add on
                    sampinds = datasample(theseinds,nOverSamp);
                    useinds=[useinds;theseinds;sampinds];
                end
                 if length(useinds)~=nMax*length(un)
                    error('mistake in resampling')
                end
            end

            % now do the classification, with the re-balanced set.           
            trnuse = trndat_thisloop(useinds,:);
            trnlabsuse = trnlabs_thisloop(useinds,:);

              %use it to predict on test set
            switch classstr
%                 case 'svmtrain_lin'
%                     obj=svmtrain(trnlabsuse,trnuse,'-t 0 -q');
%                     thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%                 case 'svmtrain_poly'
%                     obj=svmtrain(trnlabsuse,trnuse,'-t 1 -q');
%                     thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%                 case 'svmtrain_RBF'
%                     obj=svmtrain(trnlabsuse,trnuse,'-t 2 -q');
%                     thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%                 case 'svmtrain_sig'
%                     obj=svmtrain(trnlabsuse,trnuse,'-t 3 -q');
%                     thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%                 case 'fitcdiscr'
%                     obj=fitcdiscr(trnuse,trnlabsuse);
%                     thesepredlabs=predict(obj,tstdat_thisloop);
%                 case 'classify_diaglinear'
%                     thesepredlabs = classify(tstdat_thisloop,trnuse,trnlabsuse,'diagLinear');
%                 case 'classify_mahal'
%                     thesepredlabs = classify(tstdat_thisloop,trnuse,trnlabsuse,'mahalanobis');
                case 'eucDist'
                    [thesepredlabs,distAll] = eucDistClass(trnuse,tstdat_thisloop,trnlabsuse);
%                 case 'ideal_observer'
%                     [thesepredlabs,~] = ideal_observer(trnuse,tstdat_thisloop,trnlabsuse);
                case 'normEucDist'
                    [thesepredlabs,distAll] = normEucDistClass(trnuse,tstdat_thisloop,trnlabsuse);
                otherwise
                    error('your classstr does not match one of the preset options')
            end

            labs_tmp(:,bb) = thesepredlabs;
            dist_tmp(:,:,bb)  = distAll;
        end
        allpredlabs(tstinds,:) = labs_tmp;
        alldist(tstinds,:,:) = dist_tmp;
    else

        % Do classification, with the full training set.           
        trnuse = trndat_thisloop;
        trnlabsuse = trnlabs_thisloop;
        
        %use it to predict on test set
        switch classstr
%             case 'svmtrain_lin'
%                 obj=svmtrain(trnlabsuse,trnuse,'-t 0 -q');
%                 thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%             case 'svmtrain_poly'
%                 obj=svmtrain(trnlabsuse,trnuse,'-t 1 -q');
%                 thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%             case 'svmtrain_RBF'
%                 obj=svmtrain(trnlabsuse,trnuse,'-t 2 -q');
%                 thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%             case 'svmtrain_sig'
%                 obj=svmtrain(trnlabsuse,trnuse,'-t 3 -q');
%                 thesepredlabs=svmpredict(tstlabs_thisloop,tstdat_thisloop,obj);
%             case 'fitcdiscr'
%                 obj=fitcdiscr(trnuse,trnlabsuse);
%                 thesepredlabs=predict(obj,tstdat_thisloop);
%             case 'classify_diaglinear'
%                 thesepredlabs = classify(tstdat_thisloop,trnuse,trnlabsuse,'diagLinear');
%             case 'classify_mahal'
%                 thesepredlabs = classify(tstdat_thisloop,trnuse,trnlabsuse,'mahalanobis');
            case 'eucDist'
                [thesepredlabs,distAll] = eucDistClass(trnuse,tstdat_thisloop,trnlabsuse);
%             case 'ideal_observer'
%                 [thesepredlabs,~] = ideal_observer(trnuse,tstdat_thisloop,trnlabsuse);
            case 'normEucDist'
                [thesepredlabs,distAll] = normEucDistClass(trnuse,tstdat_thisloop,trnlabsuse);
            otherwise
                error('your classstr does not match one of the preset options')
        end

        % if there were more balance iterations requested, we can just
        % repeat this set of predicted labels.
        allpredlabs(tstinds,:) = repmat(thesepredlabs, 1, size(allpredlabs,2));
        alldist(tstinds,:,:) = repmat(distAll,1,1,size(alldist,3));
    end       
end

%% to get d' and accuracy: reassemble the predicted labels from all cross-validations, compute d' for entire test set at once

accAll = nan(nBalanceIter,1);
dAll = nan(nBalanceIter,1);

for bb=1:nBalanceIter

    thesepredlabs = allpredlabs(:,bb);
    if any(isnan(thesepredlabs))
        error('mistake in prediction labels')
    end
   
    accAll(bb) = mean(thesepredlabs==tstlabs);
    dAll(bb) = get_dprime(thesepredlabs,tstlabs,un);
    
end

%% get mean over the balancing iterations

acc = mean(accAll);
dprime = mean(dAll);

% get mean confidence over all balance iterations
% will be nTrials x nClasses
distEachClass = mean(alldist,3);

%% get predicted labels
% the category that the classifier assigned on most iterations

[predLabs,f,c] = mode(allpredlabs,2);

if any(cellfun('length',c)>1)
    %there was a tie, resolve randomly
    tieinds = find(cellfun('length',c)>1);
    for tt=1:length(tieinds)
        predLabs(tieinds(tt)) = datasample(c{tieinds(tt)},1);
    end
end

end