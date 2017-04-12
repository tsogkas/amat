function c = circle(r)
% CIRCLE Returnns a logical mask of size (2*r+1)x(2*r+1), in the shape of a
%   circle of radius r.
%
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: March 2017

% We do not simply use a strict equality (x^2 + y^2 == r^2) to
% avoid "gaps" beween circles of consecutive radii.
% In other words, summing circles of radii r and r+1 should
% result in a ring without holes.
r = double(r); % make sure r can take negative values
[x,y] = meshgrid(-r:r, -r:r);
c = (x.^2 + y.^2 <= r^2) & (x.^2 + y.^2 > (r-1)^2);

% Alternative definition. We prefer the previous one because it doesn't
% leave "gaps" beween circles of consecutive radii. In other words, summing
% circles of consecutive radii results in a disk mask without holes.
% c = bwperim(disk(radius));
