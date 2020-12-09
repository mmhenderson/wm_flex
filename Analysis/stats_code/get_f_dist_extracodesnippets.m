function [p_vals, ranova_table, iter] = get_f_dist(in_data, num_shuffles, make_hist)
% Makes a null distribution of f-values for the main effect(s) in a within
% subject repeated measures 1 or 2 way ANOVA. In case of a 2-way ANOVA it
% will also return an f-distribution for the interaction. Originally 
% written for a situation where I had 7 levels in one factor (ROI) and 3 in 
% the second factor (condition).
%
% OUTPUT
% p_vals: p-values associated with the real compared to shuffled f-values.
% ranova_table: the ANOVA table for the actual data (intact labels). This
% is where you can look up your data's actual F-values, and df's
% iter: the actual f-distributions in case you wanna see their shape.
%
% INPUT
% in_data: must be a matrix of num_measurements (i.e. subjects) x factor 1
% (roi) x factor 2 (condition). For a one-way ANOVA this input would 
% naturally have one dimension fewer (so only num_measurements x factor 1)
% num_shuffles: number of shuffles you wanna do (default is 1000)
% make_hist: makes histogram of your f-distributions if true (default is true)
% Written originally by Rosanne Rademaker, modified by MMH

%%
if nargin < 3, make_hist = 1; end % default is to show histogram of f's
if nargin < 2, num_shuffles = 1000; end % default is 1000 shuffles

% is this a one-way or two-way anova?
nway_indata = ndims(in_data)-1;  
assert(nway_indata==2); % STOP: can't do a one-way anova with this code.

% get some info from my input data
num_subs = size(in_data,1); % 1st dimension of in_data is number of subjects
num_f1 = size(in_data,2); % 2nd dimension is how many levels there are of factor 1
factor1_names = strseq('f1_',1:num_f1); % name the levels of my factor
num_all = num_f1; % number of all data columns combined
if nway_indata == 2 % if a two-way, also specify my second factor
    num_f2 = size(in_data,3); % 3rd dimension is how many levels there are for factor 2
    factor2_names = strseq('f2_',1:num_f2); % name the levels of my factor
    num_all = num_f1*num_f2; % number of all data columns combined
end

% make labels & within subject factor table
if nway_indata == 1 % if a one-way
    % make within subject factor table
    within_table = table(factor1_names,'VariableNames',{'factor1_names'});
elseif nway_indata == 2 % if a two-way
    % make labels for my two independent repeated measures factors
    factor1_labels = reshape(repmat(factor1_names',num_f2,1),num_f1*num_f2,1);
    factor2_labels = reshape(repmat(factor2_names',1,num_f1),num_f1*num_f2,1);
    f1_list = repelem(1:num_f1,num_f2);
    f2_list = repmat(1:num_f2,1,num_f1);
    % make within subject factor table
    within_table = table(factor1_labels, factor2_labels,'VariableNames',{'factor1_labels','factor2_labels'});
end

% make dependent variable table 
if nway_indata == 1
    dep_var = in_data;
elseif nway_indata == 2
    dep_var = reshape(permute(in_data,[1,3,2]),num_subs,num_f1*num_f2);
end
my_str = strseq('V',11:10+num_all)';
dep_var_table = array2table(dep_var,'VariableNames',my_str);

% build the model
rm = fitrm(dep_var_table,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);

% do the anova on the in_data
if nway_indata == 1
    [ranova_table] = ranova(rm,'WithinModel','factor1_names');
elseif nway_indata == 2
    [ranova_table] = ranova(rm,'WithinModel','factor1_labels*factor2_labels');
end
stats_real = ranova_table{[3,5,7],[4,6]};

%% do the shuffle
% iter = NaN(num_shuffles,3);

% stats_shuff1 = nan(num_shuffles, 3, 2);
% stats_shuff2 = nan(num_shuffles, 3, 2);
% stats_shuffint1 = nan(num_shuffles, 3, 2);
% stats_shuffint2 = nan(num_shuffles, 3, 2);
stats_shuffall = nan(num_shuffles, 3, 2);

% barcolors = plasma(4);
% barcolors = repmat(barcolors(1:2,:),14,1);
% figure;hold all;
% b=bar(mean(dep_var,1),'FaceColor','flat','CData',barcolors);
% title('real')
%%
for ii = 1:num_shuffles
        
%     dep_var_shuff_f1 = nan(size(dep_var));  % data w f1 labels shuffled
%     dep_var_shuff_f2 = nan(size(dep_var));  % data w f2 labels shuffled
    dep_var_shuff_all = nan(size(dep_var)); 
%     dep_var_shuff_f1_keepf2 = nan(size(dep_var)); 
%     dep_var_shuff_f2_keepf1 = nan(size(dep_var)); 

    % for each subject
    for n = 1:num_subs
        
%         % shuffle f2 within levels of f1
%         for ff1 = 1:num_f1
%            inds = f1_list==ff1;
%            shuf_indx = randperm(sum(inds));
%            dat2shuff = dep_var(n,inds);
%            % since we're shuffling within levels of f1, this results in the
%            % f2 labels being shuffled.
%            dep_var_shuff_f2(n,inds) = dat2shuff(shuf_indx);
%         end
%         
%         % shuffle f1 within levels of f2
%         for ff2 = 1:num_f2
%            inds = f2_list==ff2;
%            shuf_indx = randperm(sum(inds));
%            dat2shuff = dep_var(n,inds);
%            % since we're shuffling within levels of f2, this results in the
%            % f1 labels being shuffled.
%            dep_var_shuff_f1(n,inds) = dat2shuff(shuf_indx);
%         end
        
        % shuffle all levels at once
        shuf_indx = randperm(num_all);
        dep_var_shuff_all(n,:) = dep_var(n,shuf_indx);
        
%         % shuffle f1, keeping f2 groups together for each sub
%         shuf_indx = randperm(num_f1);
%         for ff2 = 1:num_f2
%            inds = f2_list==ff2;
%            dat2shuff = dep_var(n,inds);
%            % using same ordering for each level of f2
%            dep_var_shuff_f1_keepf2(n,inds) = dat2shuff(shuf_indx);
%         end
%         
%         % shuffle f2, keeping f1 groups together for each subj
%         shuf_indx = randperm(num_f2);
%         for ff1 = 1:num_f1
%            inds = f1_list==ff1;
%            dat2shuff = dep_var(n,inds);
%            % using same ordering for each level of f2
%            dep_var_shuff_f2_keepf1(n,inds) = dat2shuff(shuf_indx);
%         end
    end
    
%     figure;hold all;
%     b=bar(mean(dep_var_shuff_f1,1),'FaceColor','flat','CData',barcolors);
%     title('shuff f1 within levels of f2')
%     
%     figure;hold all;
%     b=bar(mean(dep_var_shuff_f2,1),'FaceColor','flat','CData',barcolors);
%     title('shuff f2 within levels of f1')
%     
%     figure;hold all;
%     b=bar(mean(dep_var_shuff_all,1),'FaceColor','flat','CData',barcolors);
%     title('shuff all')
%     
%     figure;hold all;
%     b=bar(mean(dep_var_shuff_f1_keepf2,1),'FaceColor','flat','CData',barcolors);
%     title('shuff f1, keeping f2 together')
    
%     dep_var_table_shuff_f1_keepf2 = array2table(dep_var_shuff_f1_keepf2,'VariableNames',my_str);
%     rm_shuff_f1_keepf2 = fitrm(dep_var_table_shuff_f1_keepf2,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);
%     [ranova_tmp] = ranova(rm_shuff_f1_keepf2,'WithinModel','factor1_labels*factor2_labels');
%     
%     stats_shuff1(ii,:,:) = ranova_tmp{[3,5,7],[4,5]};
% %     figure;hold all;
% %     b=bar(mean(dep_var_shuff_f2_keepf1,1),'FaceColor','flat','CData',barcolors);
% %     title('shuff f2, keeping f1 together')
% 
%     dep_var_table_shuff_f2_keepf1 = array2table(dep_var_shuff_f2_keepf1,'VariableNames',my_str);
%     rm_shuff_f2_keepf1 = fitrm(dep_var_table_shuff_f2_keepf1,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);
%     [ranova_tmp] = ranova(rm_shuff_f2_keepf1,'WithinModel','factor1_labels*factor2_labels');
%     
%     stats_shuff2(ii,:,:) = ranova_tmp{[3,5,7],[4,5]};
%     % run the test for the data with f1 shuffled (this will tell
%     % significance of f1's real labels)
%     dep_var_table_shuff_f1 = array2table(dep_var_shuff_f1,'VariableNames',my_str);
%     rm_shuff_f1 = fitrm(dep_var_table_shuff_f1,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);
%     [ranova_tmp] = ranova(rm_shuff_f1,'WithinModel','factor1_labels*factor2_labels');
% %     iter(i,1) = ranova_tmp{[3],4}';
% %     iter(i,3) = ranova_tmp{[7],4}';
%     stats_shuffint1(ii,:,:) = ranova_tmp{[3,5,7],[4,5]};
% 
%     % run the test for the data with f2 shuffled (this will tell
%     % significance of f2's real labels)
%     dep_var_table_shuff_f2 = array2table(dep_var_shuff_f2,'VariableNames',my_str);
%     rm_shuff_f2 = fitrm(dep_var_table_shuff_f2,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);
%     [ranova_tmp] = ranova(rm_shuff_f2,'WithinModel','factor1_labels*factor2_labels');
% %     iter(i,2) = ranova_tmp{[5],4}';
% %     iter(i,4) = ranova_tmp{[7],4}';
%     stats_shuffint2(ii,:,:) = ranova_tmp{[3,5,7],[4,5]};
%     % run the test for the data with labs shuffled (this will tell
%     % significance of all the real labels)
    dep_var_table_shuff_all = array2table(dep_var_shuff_all,'VariableNames',my_str);
    rm_shuff_all = fitrm(dep_var_table_shuff_all,[char(my_str(1)),'-', char(my_str(end)),'~1'],'WithinDesign',within_table);
    [ranova_tmp] = ranova(rm_shuff_all,'WithinModel','factor1_labels*factor2_labels');
    stats_shuffall(ii,:,:) = ranova_tmp{[3,5,7],[4,5]};
%     iter(i,:) = ranova_tmp{[3,5,7],4}';

end
%%
% p_shuff1 = mean(squeeze(stats_shuff1(:,:,1))>=repmat(squeeze(stats_real(:,1))',num_shuffles,1),1)
% p_shuff2 = mean(squeeze(stats_shuff2(:,:,1))>=repmat(squeeze(stats_real(:,1))',num_shuffles,1),1)
% p_shuffint1 = mean(squeeze(stats_shuffint1(:,:,1))>=repmat(squeeze(stats_real(:,1))',num_shuffles,1),1)
% p_shuffint2= mean(squeeze(stats_shuffint2(:,:,1))>=repmat(squeeze(stats_real(:,1))',num_shuffles,1),1)
p_shuffall = mean(squeeze(stats_shuffall(:,:,1))>=repmat(squeeze(stats_real(:,1))',num_shuffles,1),1);

% f_shuff1 = mean(squeeze(stats_shuff1(:,:,1)))
% f_shuff2 = mean(squeeze(stats_shuff2(:,:,1)))
% f_shuffint1 = mean(squeeze(stats_shuffint1(:,:,1)))
% f_shuffint2 = mean(squeeze(stats_shuffint2(:,:,1)))
f_shuffall = mean(squeeze(stats_shuffall(:,:,1)));
%%

p_vals = p_shuffall;
% if nway_indata == 1 % if a one-way
%     p_vals = sum(iter(:,1)>ranova_table{3,4})/num_shuffles;
% elseif nway_indata == 2 % if a two-way
%     p_vals = [sum(iter(:,1)>ranova_table{3,4})/num_shuffles sum(iter(:,2)>ranova_table{5,4})/num_shuffles sum(iter(:,3)>ranova_table{7,4})/num_shuffles];
% end


%% plot if you wanna
% if make_hist
%     figure;
%     if nway_indata == 1 % if a one-way
%         data_fs = ranova_table{3,4}';
%         hold on; title(['factor 1, p = ', num2str(p_vals(1))]); xlabel('F Values'); ylabel('Frequency');
%         h = histogram(iter(:,1),0:max(iter(:,1))/50:max(iter(:,1))); plot([data_fs(1) data_fs(1)],[0 max(h.Values)+max(h.Values)/10],'r');        
%     elseif nway_indata == 2 % if a two-way        
%         data_fs = ranova_table{[3 5 7],4}';
%         subplot(1,3,1); hold on; title(['factor 1, p = ', num2str(p_vals(1))]); xlabel('F Values'); ylabel('Frequency');
%         h = histogram(iter(:,1),0:max(iter(:,1))/50:max(iter(:,1))); plot([data_fs(1) data_fs(1)],[0 max(h.Values)+max(h.Values)/10],'r');
%         subplot(1,3,2); hold on; title(['factor 2, p = ', num2str(p_vals(2))]);
%         h = histogram(iter(:,2),0:max(iter(:,2))/50:max(iter(:,2))); plot([data_fs(2) data_fs(2)],[0 max(h.Values)+max(h.Values)/10],'r');
%         subplot(1,3,3); hold on; title(['interaction, p = ', num2str(p_vals(3))]);
%         h = histogram(iter(:,3),0:max(iter(:,3))/50:max(iter(:,3)));plot([data_fs(3) data_fs(3)],[0 max(h.Values)+max(h.Values)/10],'r');
%     end
% end




