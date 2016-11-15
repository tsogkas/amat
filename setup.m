%% Setup global parameters and preprocess image
numScales = 25;
numBins   = 128;
img0 = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
% img0 = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'), [128 128], 'nearest')); 
[H,W,numChannels] = size(img0);
imgLab = rgb2lab(img0);
img = clusterImageValues(imgLab, 5); % simplify input
% img = imgLab;
% img = img0;

%% Construct filters, calculate perimeneters and disk areas
filters = cell(numScales,1);
for r=1:numScales, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% Compute individual encodings for color and texture
% f = imageEncoding(img,filters);
f = imageEncoding(binImage(img,numBins),filters,'hist',numBins);

%% Compute decodings g and reconstruction errors at all points and scales
reconstructionError = imageError(img,f,filters);
reconstructionError0 = reshape(reconstructionError, H*W,numScales);