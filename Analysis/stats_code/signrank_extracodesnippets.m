%%
% Wilcoxon signed rank test
% for each permutation iteration, use this test to compare real data for all subj to
% shuffled data for all subj.
stat_iters_sr = nan(nROIs, nConds, nPermIter); 

for vv=1:nROIs
    for cc=1:nConds
        x = vals(:,vv,cc);
        for ii=1:nPermIter
            y = randvals(:,vv,cc,ii);
            % testing whether 1 (real) is greater than 2 (shuff)
            
%             diffxy = x(:) - y(:);           
%             epsdiff = eps(x(:)) + eps(y(:));
%             
%             % Remove missing data
%             t = isnan(diffxy);
%             diffxy(t) = [];
%             epsdiff(t) = [];
%            
%             t = (abs(diffxy) <= epsdiff);
%             diffxy(t) = [];
%             epsdiff(t) = [];
% 
%             n = length(diffxy);
%             chance_sr_stat = n*(n+1)/4;
%             
%             % Calculations for Sign Rank Test
% 
%             % Find positive differences and ranks of absolute differences
%             iPos = (diffxy>0);
%             [tie_rank, tieadj] = tiedrank(abs(diffxy),0,0,epsdiff);
% 
%             % Compute signed rank statistic (most extreme version)
%             w = sum(tie_rank(iPos));
%             w_pos = sum(tie_rank(iPos));
%             w_neg = sum((-1)*tie_rank(~iPos));
%             w_total = w_pos+w_neg;
%             w2 = (n*(n+1)/2) - w;
%             assert(w_total==(w-w2))
%             
% %             if w_total==0
% %                 pause
% %             end
%             
%             
%             [p2, h, stats] = signrank(x,y,'tail','right','method','approximate');
            
%             if stats.signedrank<chance_sr_stat
%                 assert(stats.zval<0)
%                 assert(w<=w2)
%                 assert(w_total<=0);
%             elseif stats.signedrank==chance_sr_stat
%                 assert(stats.zval==0)
%                 assert(w>=w2)
%                 assert(w_total>=0)
%             else
%                 assert(stats.zval>0);
%             end
            
%             if w_total<0
%                 assert(w<w2)
% %                 assert(stats.zval<0)
%                 assert(stats.signedrank<chance_sr_stat)
%             elseif w_total>0
%                 assert(w>w2)
% %                 assert(stats.zval>0)
%                 assert(stats.signedrank>chance_sr_stat)
%             elseif w_total==0
%                 assert(w==w2)
% %                 assert(stats.zval==0)
%                 assert(stats.signedrank==chance_sr_stat)
%             end
                   

            w = signrank_MMH(x,y);
            % this is a 1 if the statistic is below chance (e.g. real not
            % greater than shuffled)
            stat_iters_sr(vv,cc,ii) = w<0;

        end
    end
end

% final p value is the proportion of iterations where real was NOT
% significantly greater than shuffled (e.g. large p value)
% p_rs = mean(p_iters_rs>0.05, 3);
p_sr = mean(stat_iters_sr, 3);
% p_t = mean(p_iters_t>0.5, 3);
% is_sig=p_rs<0.05

%% more extra code for using sign rank to do paired t-tests (there are a few ways to do this)

%% now doing pairwise condition comparisons - paired t-test.
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 345455;
rng(rndseed,'twister')
% p_diff = nan(nROIs,1);
% % p_tt = nan(nROIs,1);
% real_greater_sr = nan(nROIs, nPermIter);
% rand_greater_sr = nan(nROIs, nPermIter);

real_sr_stat = nan(nROIs,1);
rand_sr_stat = nan(nROIs, nPermIter);

for vv=1:nROIs
    realvals = squeeze(vals(:,vv,:));
%     realdiff = diff(realvals,[],2);
    real_sr_stat(vv) = signrank_MMH(realvals(:,1),realvals(:,2));
%     realmeandiff = mean(realdiff);
    inds2swap = double(randn(nSubj,nPermIter)>0);
    inds2swap(inds2swap==0) = -1;
%     randmeandiffs = nan(nPermIter,1);
    parfor ii=1:nPermIter          
        
        randvals=realvals;
        randvals(inds2swap(:,ii)==-1,:) = randvals(inds2swap(:,ii)==-1,[2,1]);
        
        rand_sr_stat(vv,ii) = signrank_MMH(randvals(:,1),randvals(:,2));
        
%         % randomly swap which conditions are which, keeping subject labels
%         % same (this amounts to negating some of the differences)
%         randdiff = realdiff.*inds2swap(:,ii);
%         randmeandiffs(ii) = mean(randdiff);
%         w = signrank_MMH(realdiff,randdiff);
%         % w>0 indicates real>null, w<0 indicates real<null
%         % will use this to compute two-tailed p-value.
%         % including zero here because w=0 provides evidence in favor of the
%         % null hypothesis, and we want to count those times in both
%         % p-values.
%         real_greater_sr(vv,ii) = w>=0; 
%         rand_greater_sr(vv,ii) = w<=0;
    end

%     p_diff(vv) = 2*min([mean(realmeandiff<=randmeandiffs), mean(realmeandiff>=randmeandiffs)]);
    % also doing parametric paired t-test for comparison
%     [h,p,ci,stats]=ttest(realdiff);
%     p_tt(vv) = p;
end

% p_diff_sr = 2*min([mean(real_greater_sr,2), mean(rand_greater_sr,2)],[],2);
 

p_diff_sr = 2*min([mean(repmat(real_sr_stat,1,nPermIter)>=rand_sr_stat,2), ...
    mean(repmat(real_sr_stat,1,nPermIter)<=rand_sr_stat,2)],[],2);
p_diff = p_diff_sr;
diff_is_sig = p_diff<alpha;

% print out which areas show a significant condition effect across all subj
array2table([diff_is_sig(vis_inds), p_diff(vis_inds)],'RowNames',vis_names,'VariableNames',{'cond_diff','p'})
