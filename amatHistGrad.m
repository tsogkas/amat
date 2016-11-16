%% Setup global parameters and preprocess image
numScales = 25;
numBins   = 64;
% img0 = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
img0 = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/66075.jpg'), [128 128], 'nearest')); 
[H,W,numChannels] = size(img0);
imgLab = rgb2lab(img0);
imgClustered = clusterImageValues(imgLab, 5); % simplify input
% img = imgLab;
img = img0;

%% Construct filters, calculate perimeneters and disk areas
filters = cell(numScales,1); for r=1:numScales, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% Compute individual encodings for color and texture
f = cat(3, binImage(imgLab,numBins),textonMap(img0, numBins));
h = zeros(H,W,size(f,3),numBins,numScales);
for c=1:size(f,3)
    imgc = f(:,:,c);
    for b=1:numBins
        imgcb = double(imgc == b);
        for s=1:numScales
            D = double(filters{s})/nnz(filters{s});
            h(:,:,c,b,s) = conv2(imgcb, D, 'same');
        end
    end
end

%% Compute chi^2 distances
histGradChi = zeros(H,W,size(f,3),numScales);
for r=1:numScales-1
    histGradChi(:,:,:,r) = 0.5*sum(((h(:,:,:,:,r+1)-h(:,:,:,:,r)).^2) ./ ...
                           (h(:,:,:,:,r+1)+h(:,:,:,:,r)+eps), 4);
end

%% Compute histogram intersection 
histGradIntersect = zeros(H,W,size(f,3),numScales);
for r=1:numScales-1
    histGradIntersect(:,:,:,r) = 1-sum(min(h(:,:,:,:,r+1),h(:,:,:,:,r)),4); 
end

%% Compute max bin distance
histGradMaxBin = zeros(H,W,size(f,3),numScales);
for r=1:numScales-1
    histGradMaxBin(:,:,:,r) = max(abs(h(:,:,:,:,r+1)-h(:,:,:,:,r)),[],4);
end
