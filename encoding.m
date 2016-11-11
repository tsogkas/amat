function f = encoding(img,filters,method,numBins)
% ENCODING Compute encodings f across the image for all filters.
% 
%   f = ENCODING(img,filters,method) computes the encodings of
%   the filters contained in the cell array "filters", at all locations
%   in the image. Filters are typically disk shaped with varying radii.
%   "Method" controls the type of the encoding used, which is one of the 
%   following supported options (the one enclosed in brackets is the 
%   default method):
% 
%   {'average'}: computes the simple average of the disk-shaped region of
%                radius r, centered at the point (x,y).
%   'hist-color': computes a histogram of the binned color values in the disk
%   'hist-text' : computes a histogram of the binned texton values in the disk
%   'hist'      : computes a histogram of the binned color and texton
%                 values in the disk.
% 
%   If the chosen method is histogram-based, then an additional parameter
%   numBins {32} is required, determining the number of bins used. numBins
%   can also be a Cx1 vector, where C is the number of image channels, to
%   allow for a different number of bins per channel.
% 
%   See also: conv2
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, method  = 'average'; end
if nargin < 4, numBins = 32; end

switch method
    case 'average'
        f = meanEncodings(img,filters);
    case 'hist-color'
        img = labBins(img,numBins);
        f   = histogram(img,filters,numBins);
    case 'hist-text'
        img = textons(img,numBins);
        f   = histogram(img,filters,numBins);
    case 'hist'       
        lab   = labBins(img,numBins);
        text  = textons(img,numBins);
        flab  = histogram(lab, filters,numBins);
        ftext = histogram(text,filters,numBins);
        f     = cat(4,flab,ftext);
    otherwise, error('Method is not supported')
end


% -------------------------------------------------------------------------
function f = meanEncodings(img,filters)
% -------------------------------------------------------------------------
[H,W,numChannels] = size(img);
numScales = numel(filters);
f = zeros(H,W,numChannels,numScales);
for c=1:numChannels
    for s=1:numScales
        f(:,:,c,s) = conv2(img(:,:,c), ...
            double(filters{s})/nnz(filters{s}), 'same');
    end
end

% -------------------------------------------------------------------------
function f = histogram(img,filters,numBins)
% -------------------------------------------------------------------------
[H,W,numChannels] = size(img);
numScales = numel(filters);
f = zeros(H,W,numChannels,numBins,numScales);
for c=1:numChannels
    imgc = img(:,:,c);
    parfor b=1:numBins
        imgcb = double(imgc == b);
        for s=1:numScales
            f(:,:,c,b,s) = conv2(imgcb, ...
                double(filters{s})/nnz(filters{s}), 'same');
        end
    end
end

% -------------------------------------------------------------------------
function lab = labBins(img,numBins)
% -------------------------------------------------------------------------
if isscalar(numBins), numBins = repmat(numBins, [3,1]); end
if ismatrix(img) % grayscale image
    lab = max(1,ceil(img*numBins(1)));
else             % rgb image
    lab = rgb2lab(img);
    for i=1:3
        lab(:,:,i) = max(1,ceil(lab(:,:,i)*numBins(i)));
    end
end

% -------------------------------------------------------------------------
% The following code was taken from the Berkeley pb boundary detector package.
% https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/
% We include everything in a single file for convenience.
% -------------------------------------------------------------------------
function text = textons(img,numBins)
% -------------------------------------------------------------------------
no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
fname = sprintf('unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,numBins);
textonData = load(fname); % defines fb,tex,tsim
if ismatrix(img), tmapim = im; else tmapim = rgb2gray(img); end
text = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);

function [map] = assignTextons(fim,textons)

% Berkeley code
d = numel(fim);
n = numel(fim{1});
data = zeros(d,n);
for i = 1:d,
  data(i,:) = fim{i}(:)';
end
d2 = distSqr(data,textons);
[y,map] = min(d2,[],2);
[w,h] = size(fim{1});
map = reshape(map,w,h);

function [fim] = fbRun(fb,im)
% function [fim] = fbRun(fb,im)
%
% Run a filterbank on an image with reflected boundary conditions.
%
% See also fbCreate,padReflect.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

% find the max filter size
maxsz = max(size(fb{1}));
for i = 1:numel(fb),
  maxsz = max(maxsz,max(size(fb{i})));
end

% pad the image 
r = floor(maxsz/2);
impad = padReflect(im,r);
%  Berkeley code
% run the filterbank on the padded image, and crop the result back
% to the original image size
fim = cell(size(fb));
for i = 1:numel(fb),
  if size(fb{i},1)<50,
    fim{i} = conv2(impad,fb{i},'same');
  else
    fim{i} = fftconv2(impad,fb{i});
  end
  fim{i} = fim{i}(r+1:end-r,r+1:end-r);
end


