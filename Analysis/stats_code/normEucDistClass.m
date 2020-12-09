function [label,normEucDistAll] = normEucDistClass( train, test, group )
% Calculate the NORMALIZED Euclidean distance for each trial in "test" to each of the
% groups defined in "train"

% IN: 
    % train is [nTrialsTraining x nVoxels];
    % test is [nTrialsTesting x nVoxels];
    % group is [nTrialsTraining x 1]

% OUT:
    % label is [nTrialsTesting x 1] (this is the predicted label, as an
        % index into unique(group))
    % normEucDistAll is [nTrialsTesting x nCond] where 
        % ncond = numel(unique(group))
    
% MMH 2/13/18
%%
nconds=length(unique(group));
conds=unique(group);
nvox = size(test,2);

% first, go through each condition, get its mean, variance and number of
% samples in training set
meanrespeach = zeros(nconds,nvox);
vareach = zeros(nconds,nvox);
neach = zeros(nconds,1);
for cc=1:nconds
    % find the trials of interest in training set    
    meanrespeach(cc,:) = mean(train(group==conds(cc),:),1);
    vareach(cc,:) = var(train(group==conds(cc),:),[],1);
    neach(cc) = sum(group==conds(cc));
end
    
% use this to get the pooled variance for each voxel
pooledvar = sum((vareach.*repmat(neach-1,1,nvox)),1)/sum(neach-1);

% now loop through test set trials, and find their normalized euclidean
% distance to each of the training set conditions
normEucDistAll = zeros(size(test,1),nconds);
for cc=1:nconds
%     for tt=1:size(test,1)
        sumofsq = sum(((test-repmat(meanrespeach(cc,:),size(test,1),1))./repmat(pooledvar,size(test,1),1)).^2,2);
        normEucDistAll(:,cc) = sqrt(sumofsq);
%     end
end

% finally, assign a label to each of your testing set trials, choosing from
% the original labels in group
[~,colind]=min(normEucDistAll,[],2);
label = zeros(size(colind));
for ii=1:nconds
    label(colind==ii) = conds(ii);
end
end

