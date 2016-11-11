function f = computeDiskEncodings(img,filters,method)
% COMPUTEDISKENCODINGS Compute encodings f across the image for all filters
% 
%   f = computeDiskEncodings(img,filters,method) computes the encodings of
%   the disk-shaped filters contained in the cell array "filters", at all 
%   locations in the image. "Method" controls the type of the encoding
%   used, which is one of the following supported options (the one enclosed
%   in brackets is the default method):
% 
%   {'average'}: computes the simple average of the disk-shaped region of
%                radius r, centered at the point (x,y).
%   'hist-color': computes a histogram of the binned color values in the disk
%   'hist-text' : computes a histogram of the binned texton values in the disk
%   'hist'      : computes a histogram of the binned color and texton
%                 values in the disk.
% 
%   See also: conv2
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

[H,W,numChannels] = size(img);
numScales = numel(filters);
f = zeros(H,W,numChannels,numScales);

switch method
    case 'average'
        for c=1:numChannels
            for r=1:numScales
                f(:,:,c,r) = conv2(img(:,:,c), ...
                    double(filters{r})/nnz(filters{r}), 'same');
            end
        end
    case 'hist-color'
    case 'hist-text'
    case 'hist'
    otherwise, error('Method is not supported')
end

function b = binImage(img)

function c = computeColorHistoGrams()

function t = computeTextureHistograms()


% --- Create CIE LAB-space colormap for input image (RGB or grayscale) ----
function cmap = computeColorMap(im,nbins)
% im:       original image
% nbins:    number of bin labels

if nargin<2, nbins = ones(3,1)*32; end
if numel(nbins)==1, nbins = ones(3,1)*nbins; end

if ismatrix(im), % grayscale image
    cmap = max(1,ceil(im*nbins(1)));
else % RGB image
    % convert gamma-corrected image to LAB and scale values into [0,1]
    % min and max values for a,b channels of LAB
    % used to scale values into the unit interval
    abmin = -73;
    abmax = 95;
%   gamma = 2.5; lab = RGB2Lab(im.^gamma);
    lab = applycform(im, makecform('srgb2lab'));
    lab(:,:,1) = lab(:,:,1) ./ 100;
    lab(:,:,2) = (lab(:,:,2) - abmin) ./ (abmax-abmin);
    lab(:,:,3) = (lab(:,:,3) - abmin) ./ (abmax-abmin);
    lab(:,:,2) = max(0,min(1,lab(:,:,2)));
    lab(:,:,3) = max(0,min(1,lab(:,:,3)));
        
    cmap = zeros(size(im,1),size(im,2),3);
    for i=1:3
        cmap(:,:,i) = max(1,ceil(lab(:,:,i)*nbins(i))); 
    end
end

% --- Create texton map of an image (Berkeley Pb code) --------------------
function tmap = computeTextonMap(im,k)
% im:   original image
% k:    number of texton labels

if nargin<2, k = 64; end

no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
fname = sprintf('unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
textonData = load(fname); % defines fb,tex,tsim
if size(im,3)==3, tmapim = rgb2gray(im); else tmapim = im; end
tmap = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);

