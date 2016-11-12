function f = encodeDiskPatch(patch, method, numBins)
% ENCODEDISKPATCH returns an encoding of the pixels contained in a patch.
% 
%   p = ENCODEDISKPATCH(patch,method) returns the encoding of the patch
%   according to one of the following methods (default enclosed in 
%   brackets):
% 
%   {'average'} : simple average 
%   'hist'      : histogram of the binned values in the patch
% 
%   See also: histcounts, extractDiskPatch
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 2, method  = 'average'; end
if nargin < 3, numBins = 32; end

% Force the input into a standard shape.
if isvector(patch) % grayscale image
    patch = patch(:); 
else               % RGB image
    assert(ismatrix(patch) && size(patch,2) == 3, 'The input should be Nx3')
end 

% Compute encoding
switch method
    case 'average'
        f   = mean(patch); % scalar or 1x3 vector 
    case 'hist'
        f = zeros(numBins,size(patch,2));
        for i=1:size(patch,2)
            f(:,i) = histcounts(patch(:,i),1:numBins);
        end
    otherwise, error('Method is not supported')
end