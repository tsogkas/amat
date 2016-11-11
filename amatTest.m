%% Setup global parameters
numScales = 20;
numBins   = 32;
radius = 1:numScales;

%% Preprocess image
img = imread('google.jpg');
img = imresize(img, [128 128], 'nearest');
img = im2double(img);
img0= img; % store original RGB image (checkpoint)

img = rgb2lab(img); imgLab = img; % convert to CIE Lab color space & checkpoint
img = max(1,ceil(img * numBins)); imgBinned = img;  % bin image & checkpoint
[H,W,numChannels] = size(img);
% assert(isinrange(img, [1,numBins]))

%% Make sure that circular filters cover the whole area
% Unfortunately a percentage of the pixel IS NOT covered by using circles
% of increased radii
[x,y] = meshgrid(-radius(end):radius(end),-radius(end):radius(end));
circleFilters = zeros(size(x,1),size(x,2),numScales);
for i=1:numScales
    circleFilters(:,:,i) = bwperim(x.^2 + y.^2 <= radius(i)^2);
%     figure, imshow(circleFilters(:,:,i));
end
areaMask = sum(circleFilters(:,:,1:end),3);
% figure, imshow(areaMask)

%% Construct filters
convMaps = zeros([size(img), numScales]);
circleFilters = cell(numScales,1);
circlePerimeter = zeros(numScales,1);
for i=1:numScales
    circleFilters{i} = circle(radius(i));
    circlePerimeter(i) = nnz(circleFilters{i});
end
diskArea = cumsum(circlePerimeter);

%% Perform convolutions
circleSum = zeros(H,W,numChannels,numScales);
for c=1:numChannels
    for s=1:numScales
        circleSum(:,:,c,s) = conv2(img(:,:,c), double(circleFilters{s}), 'same');
    end
end
diskSum = cumsum(circleSum,4);

%% Test region gradients for different scales
r = 10;
ringWidth = 5;
innerDisk = diskSum(:,:,:,r)/diskArea(i);
outerRing = (diskSum(:,:,:,r+ringWidth) - diskSum(:,:,:,r)) ./ (diskArea(r+ringWidth)-diskArea(r));
% regionGrad = sum((innerDisk-outerRing).^2,3);
regionGrad = max(abs(innerDisk-outerRing),[],3);
figure, imshow(regionGrad)

%%
for ringWidth = 1:20
    regionGrad = zeros(H,W,numScales);
    for i=1:numScales-ringWidth
        innerDisk = diskSum(:,:,:,i)/diskArea(i);
        outerRing = (diskSum(:,:,:,i+ringWidth) - diskSum(:,:,:,i)) ./ (diskArea(i+ringWidth)-diskArea(i));
%         regionGrad(:,:,i) = sqrt(sum((innerDisk-outerRing).^2,3));
        regionGrad(:,:,i) = max(abs(innerDisk-outerRing),[],3);
    end
    figure, montageArray(regionGrad); title(['Region gradients, outer ring width: ' num2str(ringWidth)])
end

%% TEST GREEDY APPROACH
% Cluster image values to simplify input
numClusters = 6;
img = reshape(img0, H*W,numChannels);
[clusterInd, clusterCenters] = kmeans(img,numClusters,'MaxIter',1000);
for k=1:numClusters
    img(clusterInd == k,:) = repmat(clusterCenters(k,:), [nnz(clusterInd == k),1]);
end

% Remove grey areas around boundaries
minCluster = 1; minCount = nnz(clusterInd == 1); % Find smallest cluster
[~,maxCluster] = min(sum(bsxfun(@minus,clusterCenters,[1 1 1]).^2 ,2));
for i=2:numClusters
    if nnz(clusterInd == i) < minCount
        minCluster = i;
        minCount = nnz(clusterInd == i);
    end
end
img(clusterInd == minCluster,:) = ...
    repmat(clusterCenters(maxCluster,:), [nnz(clusterInd == minCluster),1]);
assert(size(unique(img,'rows'),1) == numClusters - 1);
img = reshape(img,H,W,numChannels);
imgClustered = img; % checkpoint


% Let's see now!
% 1. Compute encodings f(D_I(x,y,r)) at every point.
img = imgClustered;
circleSum = zeros(H,W,numChannels,numScales);
for c=1:numChannels
    for s=1:numScales
        circleSum(:,:,c,s) = conv2(img(:,:,c), double(circleFilters{s}), 'same');
    end
end
f = cumsum(circleSum,4); % (x,y,k,r)
f = bsxfun(@rdivide,f,reshape(diskArea,1,1,1,[]));

% Compute decodings g (D*_I(x,y,r)) at every point, compare with respective
% image patch D_I(x,y,r) and compute their distance.
% d = zeros(H,W,numScales);
% [x,y] = meshgrid(1:H,1:W);
% img = reshape(img, H*W, []);
% for r=1:numScales
%     disp(['Scale ' num2str(r)])
%     for i=1:H
%         for j=1:W
%             xc = j; yc = i;
%             mask = (x-xc).^2 + (y-yc).^2 <= r^2;
%             D_I = img(mask,:);
%             d(i,j,r) = sum(sum(bsxfun(@minus, D_I, reshape(f(i,j,:,r), [], numChannels)).^2));
%         end
%     end
% end

% Faster  version
[x,y] = meshgrid(1:H,1:W);
reconstructionError = zeros(H*W,numScales);
for r=1:numScales
    disp(['Computing reconstruction error for scale ' num2str(r)])
    % pad img and coordinate matrices so that im2col gives correct #nhoods
    % pixels outside the original image domain are not considered
    img      = padarray(imgClustered, [r,r], 0, 'both'); 
    xpatches = padarray(x, [r,r], inf, 'both');  
    ypatches = padarray(y, [r,r], inf, 'both');
    xpatches = im2col(xpatches,   [2*r+1, 2*r+1],'sliding');
    ypatches = im2col(ypatches,   [2*r+1, 2*r+1],'sliding');
    rpatches = im2col(img(:,:,1), [2*r+1, 2*r+1],'sliding');
    gpatches = im2col(img(:,:,2), [2*r+1, 2*r+1],'sliding');
    bpatches = im2col(img(:,:,3), [2*r+1, 2*r+1],'sliding');
    assert(size(xpatches,2) == H*W)
    assert(size(ypatches,2) == H*W)
    % subtract disk center coordinates from nhood coordinates
    xpatches = bsxfun(@minus,xpatches,x(:)');
    ypatches = bsxfun(@minus,ypatches,y(:)');
    % Find which points are inside the disk for each nhood
    indisks  = xpatches.^2 + ypatches.^2 <= r^2;
    % Compute error between original and reconstructed patches
    % WARNING!: to compute the reconstruction error correctly, we must make 
    % sure that we count only pixels INSIDE the disk
    imgpatches = cat(3,rpatches,gpatches,bpatches);
    recpatches = reshape(f(:,:,:,r), 1,[],numChannels); % img patch approximation
    recerror   = bsxfun(@minus, imgpatches, recpatches);
    recerror   = bsxfun(@times, recerror, double(indisks));
    reconstructionError(:,r) = sum(sum(recerror.^2, 3), 1);   
end

