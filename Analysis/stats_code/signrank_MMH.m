function [w] = signrank_MMH(x,y)
      
% simplified version of signrank.m 
% returns only the raw test-statistic, keeping its original sign. The sign of
% this test statistic tells you which of the groups is bigger (positive
% value means x>y, negative value means x<y)
% note: to compute the significance of this statistic on its own, would need to
% normalize by its variance, which is not done in this code. 

% requires the functions eps.m and tiedrank.m (stats toolbox)
% MMH 11/30/20

    %%

    diffxy = x(:) - y(:);           
    epsdiff = eps(x(:)) + eps(y(:));
    assert(~any(isnan(diffxy)));

    % Calculations for Sign Rank Test
    % Find positive differences and ranks of absolute differences
    iPos = (diffxy>0);  
    iNeg = (diffxy<0);
    [tie_rank, ~] = tiedrank(abs(diffxy),0,0,epsdiff);

    % Compute signed rank statistic 
    % w = sum over n (sign of the difference * rank of the difference)
    w = sum(tie_rank(iPos)) - sum(tie_rank(iNeg));
      