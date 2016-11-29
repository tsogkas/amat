function d = histogramDistance(h1,h2,method,r)
% HISTOGRAMDISTANCE Compute the distance between histograms using various
% metrics (assuming that h1 and h2 have been normalized in [0,1]).
% 
%   d = HISTOGRAMDISTANCE(h1,h2,method) computes the distance of histograms
%   h1 and h2. h1 must have the same size and dimensions and they can be:
%   
%   1) vectors: representing two histograms with B = length(h1) bins.
%   2) CxB matrices: representing C histograms (e.g. one histogram per 
%      image channel) of B bins each.
%   3) HxWxCxB arrays: where h1 and h2 correspond to histograms of B bins
%      for each pixel in a HxW image with C channels.
% 
%   Method is a string that can take one of the following values: 
% 
%   chi2: chi-squared histogram distance.
%   chi2-kernel: chi-squared histogram distance + gaussian kernel density 
%       estimate. In this case the used must provide a fourth parameter which
%       is the sigma of the gaussian kernel.
%   intersection: 1 - intersection(h1,h2) (since histogram intersection 
%       models histogram *similarity*, not distance.
% 
%   See also: histogram, histcounts
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, method = 'chi2'; end

assert(isequal(size(h1),size(h2)), 'Histogram dimensions missmatch')
if isvector(h1), [~,dim] = max(size(h1)); % works for column and row vectors
elseif ndims(h1) == 5, dim = 4;  % HxWxCxBxR
else dim = ndims(h1); % in all other cases, bins are the last dim
end

switch method
    case 'intersection'
        d = 1-sum(min(h1,h2),dim); % 1-histogram intersection
    case 'chi2'
        d = 0.5*sum((h1-h2).^2 ./ (h1 + h2 + eps), dim);
    case 'chi2-kernel'
        assert(nargin == 4 && ~isempty(r), ...
            'You must provide the sigma for the gaussian kernel')
        sz = size(h1);
        B = sz(dim);
        % Compute color bin weighted distance using gaussian kernel
        binCenters = ((1:B)-0.5)/B;
        [x,y] = meshgrid(binCenters,binCenters);
        % Compute distance at a single scale r (sigma of the gaussian)
        h1 = reshape(h1,[],B);
        h2 = reshape(h2,[],B);
        binCenterDistance = 1-exp(-(x-y).^2./(2*r.^2)); % BxB
        dabs = abs(h1-h2); % H*W*C x B
        d = sum((dabs*binCenterDistance) .* dabs, 2); % H*W*C x 1
        d = reshape(d, [sz(1:dim-1) sz(dim+1:end)]);  % HxWxC
    otherwise, error('Histogram comparison method not supported')
end