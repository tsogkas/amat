function y = histinv(h,n)
% HISTINV Sample from a distribution described by a given histogram,
% 
%   y = HISTINV(h,n) takes as input the histogram distribution h and
%   returns n samples drawn from the distribution. HISTINV uses the inverse
%   transform sampling or Smirnov transform to draw the samples. h is a 1xB
%   vector, or real numbers in [0,1], representing probabilities where B is 
%   the number of bins used.
% 
%   See also: rand
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016


if isvector(h)
    B = length(h);
    binCenters = (0:(B-1))/B;
    y = rand(n,1);
    [~,y] = max(bsxfun(@gt,cumsum(h),y),[],2);
    y = binCenters(y);
else
    % We assume for convenience that the histogram 
    % dimension is the last dimension of the input array.
    dim = ndims(h);
    sz = size(h);
    B = sz(end);
    binCenters = (0:(B-1))/B;
    y = rand([sz, 1, n]);
    [~,y] = max(bsxfun(@gt,cumsum(h,dim),y),[],dim);
    y = reshape(binCenters(y), [sz(1:dim-1),b]);
end
