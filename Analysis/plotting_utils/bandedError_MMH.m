function h = bandedError_MMH(x, y,se, color, alpha)
%h = bandedError(x, y, se, plotHandle, alpha)
%js (serences@salk.edu) 07.03.2006...inspired by some data presented
%       by clay curtis
% edited by MMH so you can set a color for the error bars

if nargin<4
    color = [0,0.3, 0.6];
    alpha = 0.3;
end
nLines = size(x,1);

assert(all(size(x)==size(y)) & all(size(x)==size(se)));

if size(x,1)>size(x,2) % column vectors
    x = x';
    y = y';
    se = se';
end

for l=1:nLines

tmpx = [x(l,:),fliplr(x(l,:))];                                   %wrap around the xData for passing to patch
tmpy = [y(l,:)+se(l,:),fliplr(y(l,:)-se(l,:))];                             %make the y-data (again, 'wrapping it around')
h = patch(tmpx, tmpy, color);
set(h, 'FaceAlpha', alpha);
set(h, 'EdgeColor','none');

end