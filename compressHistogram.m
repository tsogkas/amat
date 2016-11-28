function  y = compressHistogram(hin,method,dim)
% COMPRESSHISTOGRAM Summarizes a histogram encoding using a single value 
% (or vector).
% 
%   The supported methods are the following:
% 
%   {'mode'}: Summarizes the histogram using the centroid of the bin with 
%             maximum probability.
%   'expectation':  Weighted sum of the histogram centroids.
% 
%   See also: histogram
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, dim = 1; end
if nargin < 2, method = 'mode'; end

B = size(hin,dim);
binCenters = (0:(B-1))/B;
switch method
    case 'mode'
        [~,maxInds] = max(hin,[],dim);
        y = squeeze(binCenters(maxInds));
    case 'expectation'
        % this can become more efficient with reshape and matrix
        % multiplication but the code will won't probably generalize so
        % easily
        reshapeVector = ones(size(hin)); reshapeVector(dim) = B;
        binCenters = reshape(binCenters, reshapeVector);
        y = squeeze(sum(bsxfun(@times,hin,binCenters),dim)); 
    otherwise, error('Method not supported')
end