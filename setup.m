%% Setup global parameters and preprocess image
numScales = 25;
img0 = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
% img0 = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'), [128 128], 'nearest')); 
[H,W,numChannels] = size(img0);
% imgLab = rgb2lab(img0);
% img = clusterImageValues(imgLab, 5); % simplify input
% img = imgLab;
% img = img0;

%% Construct filters, calculate perimeneters and disk areas
filters = cell(numScales,1);
for r=1:numScales, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
f = encoding(img,filters);

%% Compute decodings g and reconstruction errors at all points and scales
[x,y] = meshgrid(1:W,1:H); % H x W
reconstructionError = zeros(H*W,numScales);
for r=1:numScales
    disp(['Computing reconstruction error for scale ' num2str(r)])
    % pad img and coordinate matrices so that im2col gives correct #nhoods
    % pixels outside the original image domain are not considered
    imgp     = padarray(img, [r,r], 0, 'both'); 
    xpatches = [repmat((-r+1):0,    [H,1]), x, repmat((W+1):(W+r),    [H,1])];
    ypatches = [repmat(((-r+1):0)', [1,W]); y; repmat(((H+1):(H+r))', [1,W])];
    xpatches = padarray(xpatches, [r,0], 'replicate', 'both');
    ypatches = padarray(ypatches, [0,r], 'replicate', 'both');  
    xpatches = im2col(xpatches,   [2*r+1, 2*r+1],'sliding'); % #nhoodPixels x H*W
    ypatches = im2col(ypatches,   [2*r+1, 2*r+1],'sliding');
    rpatches = im2col(imgp(:,:,1),[2*r+1, 2*r+1],'sliding');
    gpatches = im2col(imgp(:,:,2),[2*r+1, 2*r+1],'sliding');
    bpatches = im2col(imgp(:,:,3),[2*r+1, 2*r+1],'sliding');
    assert(size(xpatches,2) == H*W)
    assert(size(ypatches,2) == H*W)
    % subtract disk center coordinates from nhood coordinates
    xpatches = bsxfun(@minus,xpatches,x(:)');
    ypatches = bsxfun(@minus,ypatches,y(:)');
    % Find which points are inside the disk for each nhood
    indisk  = double(xpatches.^2 + ypatches.^2 <= r^2); % #nhoodPixels x H*W
    % Compute error between original and reconstructed patches
    % WARNING!: to compute the reconstruction error correctly, we must make 
    % sure that we count only pixels INSIDE the disk
    imgpatches = cat(3,rpatches,gpatches,bpatches); % #nhoodPixels x H*W x 3
    recpatches = reshape(f(:,:,:,r), 1,[],numChannels);  % g(x) 1 x H*W x 3 
    nomin = bsxfun(@minus, imgpatches, recpatches); % g(x)-f(x) #nhoodPixels x H*W x 3
    nomin = bsxfun(@times, nomin, indisk); % mask out non-disk pixels
    nomin = sum(sum(nomin.^2, 3), 1); % SUM(||g(x)-f(x)||^2)
    denom = sum(sum(bsxfun(@times, imgpatches, indisk).^2, 3), 1); % 1xH*W
    reconstructionError(:,r) = sqrt(nomin ./ denom); % squared normalized rms error  
    assert(isinrange(recpatches, [0,1+1e-8]))
    assert(isinrange(imgpatches, [0,1+1e-8]))
end
assert(isinrange(reconstructionError, [0,1]))
reconstructionError0 = reconstructionError;
clear xpatches ypatches rpatches gpatches bpatches recpatches imgpatches indisk nomin denom
