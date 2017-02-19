function c = circle(r)
r = double(r); % make sure r can take negative values
[x,y] = meshgrid(-r:r, -r:r);
c = x.^2 + y.^2 <= r^2 & x.^2 + y.^2 > (r-1)^2;

% Alternative definition. We prefer the previous one because it doesn't
% leave "gaps" beween circles of consecutive radii. In other words, summing
% circles of consecutive radii results in a disk mask without holes.
% c = bwperim(disk(radius));