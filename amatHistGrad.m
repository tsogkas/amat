%% Setup global parameters and preprocess image
R = 25; % #scales
B = 64; % #bins
errorType = 'se';
encodingType = 'average';
colorWeights = [];

% imgRGB = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
imgRGB = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/66075.jpg'), [128 128], 'nearest')); 
[H,W,C] = size(imgRGB);
imgLab = rgb2labNormalized(imgRGB);
% imgClustered = clusterImageValues(imgLab, 5); % simplify input
if strcmp(errorType, 'dssim'), img = imgRGB; else img = imgLab; end

%% Construct filters, calculate perimeneters and disk areas
filters = cell(R,1); for r=1:R, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% Compute individual encodings for color and texture
f = cat(3, binImage(imgLab,B),textonMap(imgRGB, B));
h = zeros(H,W,size(f,3),B,R);
for c=1:size(f,3)
    imgc = f(:,:,c);
    for b=1:B
        imgcb = double(imgc == b);
        for s=1:R
            D = double(filters{s})/nnz(filters{s});
            h(:,:,c,b,s) = conv2(imgcb, D, 'same');
        end
    end
end

%% Compute chi^2 distances
histGradChi = zeros(H,W,size(f,3),R);
for r=1:R-1
    histGradChi(:,:,:,r) = 0.5*sum(((h(:,:,:,:,r+1)-h(:,:,:,:,r)).^2) ./ ...
                           (h(:,:,:,:,r+1)+h(:,:,:,:,r)+eps), 4);
end

%% Compute histogram intersection 
histGradIntersect = zeros(H,W,size(f,3),R);
for r=1:R-1
    histGradIntersect(:,:,:,r) = 1-sum(min(h(:,:,:,:,r+1),h(:,:,:,:,r)),4); 
end

%% Compute max bin distance
histGradMaxBin = zeros(H,W,size(f,3),R);
for r=1:R-1
    histGradMaxBin(:,:,:,r) = max(abs(h(:,:,:,:,r+1)-h(:,:,:,:,r)),[],4);
end

%% Debug
r2 = 2;
r1 = 1;
c  = 1;
figure(1);
subplot(131); imagesc(0.5*sum(((h(:,:,c,:,r2)-h(:,:,c,:,r1)).^2) ./ (h(:,:,c,:,r2)+h(:,:,c,:,r1)+eps), 4),[0,1]); title('Chi-square ')
subplot(132); imagesc(max(abs(h(:,:,c,:,r2)-h(:,:,c,:,r1)),[],4),[0,1]); title('Max absolute')
subplot(133); imagesc(1-sum(min(h(:,:,c,:,r2),h(:,:,c,:,r1)),4),[0,1]); title('Histogram intersection')