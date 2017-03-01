function out = bwthin(in,n)
if nargin < 2, n = inf; end
out = bwmorph(in,'thin',n);