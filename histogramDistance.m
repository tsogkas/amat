function d = histogramDistance(h1,h2,method,params)
% HISTOGRAMDISTANCE Compute the distance between histograms using various
% metrics (assuming that h1 and h2 have been normalized in [0,1]).
% 
%   d = HISTOGRAMDISTANCE(h1,h2,method) computes the distance of histograms
%   h1 and h2. 
%   
%   d = HISTOGRAMDISTANCE(h1,h2,method,params) also allows the user to give
%   method-specific parameters as input.
% 
%   h1 must have the same size and dimensions and they can be:
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

% TODO: maybe squeeze cases with sum() and max()?
% TODO: add chi2-gaussian as sub-case of chi-quadratic?
switch method
    case {'L0','hellinger'} % in [0,1]
        d = sqrt(sum((sqrt(h1)-sqrt(h2)).^2, dim)) / sqrt(2);
    case {'L1','manhattan'} % in [0,1]
        d = sum(abs(h1-h2),dim);
    case {'L2','euclidean'} % in [0,1]
        d = sqrt(sum((h1-h2).^2, dim));
    case {'Linf','chebyshev','KS','ks','kolmogorov-smirnov'} % in [0,1]
        d = squeeze(max(abs(h1-h2),[],dim)); 
    case {'Lp','fractional'} % in [0,1]
        p = params;
        d = sum(abs((h1-h2).^p), dim).^(1/p);
    case {'KL','kl','kullback-leibler'}
        d = sum(h1 .* log( h1 ./ (h2 + eps)), dim);
    case {'bhattacharaya','bhattacharyya'} % in [0,inf)
        d = -log(sum(sqrt(h1 .* h2), dim));
    case {'cm','CM','cramer','cramer-von-mises','cramer-vm'} % need to check this out
        % d = sum((cumsum(h1,dim) - cumsum(h2,dim)).^2, dim); 
    case {'emd','EMD','earth-mover'}
        % NOT IMPLEMENTED YET
        error('Not implemented yet');
    case 'cosine'
        d = 1-sum(h1 .* h2, dim);
    case 'canberra'
        d = sum(abs(h1-h2) ./ (h1+h2+eps), dim);
    case 'pearson'
        % NOT IMPLEMENTED YET
        error('Not implemented yet');
    case 'square-chord'
        d = sum((sqrt(h1)-sqrt(h2)).^2, dim);
    case 'jeffrey' % (symmetric KL-div)
        d = sum(h1 .* log( h1 ./ (h2 + eps)) + h2 .* log( h2 ./ (h1 + eps)), dim);
    case 'jensen-shannon' 
        d = 0.5*sum(h1 .* log( 2*h1 ./ (h1+h2+eps)) + h2 .* log( 2*h2 ./ (h1+h2+eps)), dim);
    case 'quadratic'
        assert(nargin == 4 && ~isempty(params), ...
            ['You must provide the bin similarity matrix or the ' ...
             'function handle (kernel) to be applied on the bin distances'])        
        sz = size(h1);
        B = sz(dim);
        if ismatrix(params), similarityMatrix = params;
        elseif isa(params,'function_handle')
            kernel = params;
            binCenters = ((1:B)-0.5)/B;
            [x,y] = meshgrid(binCenters,binCenters);
            similarityMatrix = kernel(x,y);
        else error('You must provide either a similarity matrix or a function handle')
        end
        h1 = reshape(h1,[],B);
        h2 = reshape(h2,[],B);
        df = h1-h2;
        d  = sqrt(sum((df * similarityMatrix) .* df, 2));
        d  = reshape(d, [sz(1:dim-1) sz(dim+1:end)]); 
    case 'chi-quadratic'
        % NOT IMPLEMENTED YET
        error('Not implemented yet');
    case 'intersection'
        d = 1-sum(min(h1,h2),dim); % 1-histogram intersection
    case 'chi2'
        d = 0.5*sum( ((h1-h2).^2) ./ (h1 + h2 + eps), dim);
    case 'chi2-gaussian'
        assert(nargin == 4 && ~isempty(params), ...
            'You must provide the sigma for the gaussian kernel')
        r = params;
        sz = size(h1);
        B = sz(dim);
        % Compute color bin weighted distance using gaussian kernel
        binCenters = ((1:B)-0.5)/B;
        [x,y] = meshgrid(binCenters,binCenters);
        % Compute distance at a single scale r (sigma of the gaussian)
        h1 = reshape(h1,[],B);
        h2 = reshape(h2,[],B);
        binCenterDistance = 1-exp(-((x-y).^2) ./ (2*r^2)); % BxB
        dabs = abs(h1-h2); % H*W*C x B
        d = sqrt(sum((dabs*binCenterDistance) .* dabs, 2)); % H*W*C x 1
        d = reshape(d, [sz(1:dim-1) sz(dim+1:end)]);  % HxWxC
    otherwise, error('Histogram comparison method not supported')
end