function [label,distAll] = eucDistClass( train, test, group )
% Use a Euclidean distance classifer to predict the labels of test data
% Also output the Euclidean distance for each sample to each group

% inputs:
    % train is [ntrials_trn x nvox]
    % test is [ntrials_tst x nvox]
    % group is [ntrials_trn x 1]
    
% outputs
    % label is [ntrials_tst x 1]
        % the values in label correspond to the original values in group
    % distAll is [ntrials_tst x numel(unique(group))]
        % describes the euclidiean distance from the test trial pattern to 
        % the mean for each condition - each column corresponds to an 
        % element in group, listed in ascending sorted order.
        
% MMH 1/16/2018

%%
nconds=length(unique(group));
conds=unique(group);
% npred=size(train,2);

distAll = zeros(size(test,1),nconds);

for cc=1:nconds
    %compute average resp to this cond
    thisaverage = mean(train(group==conds(cc),:),1);

    % subtract average pattern from each test trial, square, sum, sqrt
    distAll(:,cc) = sqrt(sum((test-repmat(thisaverage,size(test,1),1)).^2,2));
   
end


[~,colind]=min(distAll,[],2);
label = zeros(size(colind));
for ii=1:nconds
    label(colind==ii) = conds(ii);
end
end

