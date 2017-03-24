function d = disk(r)
% DISK Returns a logical mask of size (2*r+1)x(2*r+1), in the shape of a
%   disk of radius r.
%
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: March 2017

r = double(r); % make sure r can take negative values
[x,y] = meshgrid(-r:r, -r:r);
d = x.^2 + y.^2 <= r^2;
