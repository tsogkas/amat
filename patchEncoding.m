function f = patchEncoding(patch, method, B)
% PATCHENCODING returns an encoding of the pixels contained in a patch.
% 
%   p = PATCHENCODING(patch,method) returns the encoding of the patch
%   according to one of the following methods (default enclosed in 
%   brackets):
% 
%   {'average'} : simple average (returns a 1xC) vector.
%   'hist'      : histogram of the binned values in the patch (returns a
%                 CxB matrix of histograms for each channel, where B is the
%                 number of bins used and is specified by the user as a 3rd
%                 input argument).
% 
%   Patch can be a Nx1 vector where N is the number of the pixels in the 
%   patch (grayscale image), a NxC matrix, where C is the number of the 
%   image channels, or a HxWxC array, where HxW are the dimensions of the 
%   patch.
%   
%   See also: histcounts, diskPatch, imageEncoding
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 2, method  = 'average'; end
if nargin < 3, B = 32; end

% Force the input into a standard shape.
if isvector(patch) % grayscale image
    patch = patch(:);
elseif ismatrix(patch) % do nothing
elseif ndims(patch) == 3
    patch = reshape(patch,[],size(patch,3));
else error('Input should be a vector, a matrix, or a 3D array.')
end 

% Compute encoding
switch method
    case 'average'
        f = mean(patch); % scalar or 1x3 vector 
    case 'hist'
        [N,C] = size(patch);
        f = zeros(C,B);
        binEdges = (0:B)/B;
        for i=1:C
            f(i,:) = histcounts(patch(:,i),binEdges);
        end
        f = f/N;
    otherwise, error('Method is not supported')
end