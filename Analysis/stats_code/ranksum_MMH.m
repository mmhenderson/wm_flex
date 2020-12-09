function [w] = ranksum_MMH(x,y)
      
% simplified version of ranksum.m 
% returns only the raw test-statistic. 
% note: to compute the significance of this statistic on its own, need to
% do more calculations not done in this code.

% requires the functions tiedrank.m (stats toolbox)
% MMH 11/30/20

%%
    % Remove missing data
    x = x( ~isnan(x) );
    y = y( ~isnan(y) );

    nx = numel(x);
    
    % Calculations for Rank Sum Test
    x = x(:);   % ensure columns
    y = y(:);

    % Compute the rank sum statistic based on the first sample
    [ranks, ~] = tiedrank([x;y]);
    srank = ranks(1:nx);
    w = sum(srank);
      