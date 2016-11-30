function f = imageEncoding(img,filters,method,B)
% IMAGEENCODING Compute encodings f across the image for all given filters. 
%   Given an input image, the encoding f(x,y,s) is a "summary" of the 
%   appearance content of the image at a region filters{s}, centered at 
%   pixel (x,y), at scale r.
% 
%   f = IMAGEENCODING(img,filters,method) computes the encodings of
%   the filters contained in the cell array "filters", at all locations
%   in the image. Filters are typically disk shaped with varying radii.
%   "Method" controls the type of the encoding used, which is one of the 
%   following supported options (default enclosed in brackets):
% 
%   {'average'} : computes the simple average of the disk-shaped region of
%                radius r, centered at the point (x,y).
%   'hist'      : computes a histogram of the binned values in the disk.
% 
%   If the chosen method is histogram-based, then an additional parameter
%   B {32} is required, determining the number of bins used. numBins
%   can also be a Cx1 vector, where C is the number of image channels, to
%   allow for a different number of bins per channel.
% 
%   NOTE: we treat channels a* and b* in the Lab color space as independent
%   for simplicity (in reality they are not).
% 
%   See also: patchEncoding
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, method  = 'average'; end
if nargin < 4, B = 32; end

switch method
    case 'average'
        f = meanEncoding(img,filters); % HxWxCxR
    case 'hist'
        f = histogramEncoding(img,filters,B); % HxWxCxBxR
    otherwise, error('Method is not supported')
end


% -------------------------------------------------------------------------
function f = meanEncoding(img,filters)
% -------------------------------------------------------------------------
[H,W,C] = size(img);
R = numel(filters);
f = zeros(H,W,C,R);
for c=1:C
    for s=1:R
        D = double(filters{s})/nnz(filters{s});
        f(:,:,c,s) = conv2(img(:,:,c), D, 'same');
    end
end

% -------------------------------------------------------------------------
function f = histogramEncoding(img,filters,B)
% -------------------------------------------------------------------------
[H,W,C] = size(img);
R = numel(filters);
f = zeros(H,W,C,B,R);
for s=1:R
    D = double(filters{s})/nnz(filters{s});
    for b=1:B
        imgb = double(img == b);
        for c=1:C
            f(:,:,c,b,s) = conv2(imgb(:,:,c), D, 'same');
        end
    end
    csum = sum(f(:,:,:,2:end,s),4);
    f(1:s,:,:,1,s) = 1-csum(1:s,:,:);
    f(:,1:s,:,1,s) = 1-csum(:,1:s,:);
    f(end-s+1:end,:,:,1,s) = 1-csum(end-s+1:end,:,:);
    f(:,end-s+1:end,:,1,s) = 1-csum(:,end-s+1:end,:);
end
assert(allvec(sum(f,4)==1),'Histograms do not sum up to 1')






