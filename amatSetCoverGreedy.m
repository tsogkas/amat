%% Setup global parameters and preprocess image
numScales = 25;
numBins   = 64;
errorType = 'se';
encodingType = 'average';
colorWeights = [.7 .15 .15];

% img0 = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
imgRGB = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/66075.jpg'), [128 128], 'nearest')); 
[H,W,numChannels] = size(imgRGB);
imgLab = rgb2labNormalized(imgRGB);
% imgClustered = clusterImageValues(imgLab, 5); % simplify input
if strcmp(errorType, 'dssim'), img = imgRGB; else img = imgLab; end

%% Construct filters, calculate perimeneters and disk areas
filters = cell(numScales,1); for r=1:numScales, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
if strcmp(encodingType,'average')
    f = imageEncoding(img,filters);
elseif strcmp(encodingType,'hist')
    f = imageEncoding(binImage(img,numBins),filters,'hist',numBins);
end

%% Compute decodings g and reconstruction errors at all points and scales
reconstructionError = imageError(img,f,filters,errorType,colorWeights);

%% Greedy approximation of the weighted set cover problem associated with AMAT
% TODO:
% - Manually compare errors and find out why after some point excessively
%   large disks are preferred than smaller disks with zero-error. 
% - Keep track of disk inclusions using sparse matrix and replace for loops
% - Maybe add an extra regularization term that discourages disks with
%   radius that does not agree with the radii of neighboring/enclosed disks

%% Initializations
amat.input          = imgRGB;
amat.reconstruction = zeros(H*W,numChannels);
amat.axis           = zeros(H,W,numChannels);
amat.radius         = zeros(H,W);
amat.depth          = zeros(H,W); 
amat.covered        = false(H,W);
% Used in the greedy approximate algorithm. Not sure how we will exploit it
amat.price          = inf(H,W); 

% Easy way to compute the number of NEW pixels that will be covered by each 
% disk if it is added in the solution, taking into account the fact that
% larger disks exceed image boundaries.
numNewPixelsCovered = ones(H,W,numScales);
for r=1:numScales
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), double(filters{r}),'same');
end
numNewPixelsCovered = reshape(numNewPixelsCovered ,H*W,numScales);

%% Error balancing and visualizaion of top (low-cost) disks
% We must *add* a scale-based regularization term, to favour larger radii
% even when the errors are 0. Dividing by the respective radius would not
% work in that case.
diskCost = reshape(reconstructionError, H*W,numScales);
% diskCost = bsxfun(@plus, diskCost, (0)./(1:numScales));
diskCostEffective = diskCost./ numNewPixelsCovered;
% Sort costs in ascending order and visualize top disks.
[sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
top = 1e2;
[yy,xx,rr] = ind2sub([H,W,numScales], indSorted(1:top));
figure(1); imshow(amat.input); 
viscircles([xx,yy],rr,'Color','w','LineWidth',0.5);
viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
title(sprintf('W: Top-%d disks, B: Top-1 disk',top))


%% Run the greedy algorithm
[x,y] = meshgrid(1:W,1:H);
while ~all(amat.covered(:))
    % Find the most cost-effective set in the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    % Build set D on the fly
    [yc,xc,rc] = ind2sub([H,W,numScales], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    newPixelsCovered = D & ~amat.covered;
    
    % Update AMAT
    amat.price(newPixelsCovered) = minCost; % Set price for new elements
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = repmat(f(yc,xc,:,rc), [nnz(newPixelsCovered),1]);
    amat.axis(yc,xc,:) = f(yc,xc,:,rc);
    amat.radius(yc,xc) = rc;

    % TL;DR: Subtract newly covered pixels from all overlapping disks. 
    % LONG STORY: We have decided to place a disk centered at (xc,yc)
    % with radius r. This disks covers all points (xd,yd) that lie within
    % radius r. For each point (xd,yd) we must find how many disks it is
    % covered by and subtract 1 from the number of NEW pixels covered by
    % those disks (if they were to be selected as part of the solution).
    % Each (xd,yd) is covered by all disks that have centers at distance <=
    % maxRadius from them.
    % find coordinates of all newly covered pixels inside D
    [yd,xd] = ind2sub([H,W],find(newPixelsCovered)); 
    for i=1:numel(yd)
        distFromCoveredPointSquared = (x-xd(i)).^2 + (y-yd(i)).^2;
        % All the points that lie within max radius distance from the 
        % current point (xd,yd) are centers of (at least) one disk that contains it.
        centersOfCoveringDisks = distFromCoveredPointSquared <= numScales^2;
        % All disks of radius >= minCoveringDistance cover (xd,yd) so we 
        % have to substract 1 from the the number of NEW  pixels that will
        % be covered by each disk if it's added in the solution.
%         minCoveringDistanceSquared = ceil(distFromCoveredPointSquared(centersOfCoveringDisks));
        for r=1:numScales
            inds = centersOfCoveringDisks & (distFromCoveredPointSquared <= r^2);
            numNewPixelsCovered(inds, r) = numNewPixelsCovered(inds, r) - 1;
            diskCost(inds,r) = diskCost(inds,r) - sum((amat.input(yd(i),xd(i),:) - f(yd(i),xd(i),:,r)).^2);
        end        
    end
    
    % Update cost effectiveness score
    diskCostEffective = diskCost./ numNewPixelsCovered;
    assert(all(isinf(diskCostEffective(sub2ind([H,W],yc,xc), 1:rc))))

    
    % Visualize progress
    if 1
        % Sort costs in ascending order to visualize updated top disks.
        [sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
        [yy,xx,rr] = ind2sub([H,W,numScales], indSorted(1:top));
        subplot(221); imshow(amat.input); 
        viscircles([xc,yc],rc, 'Color','k','EnhanceVisibility',false);
        title('Selected disk');
        subplot(222); imshow(bsxfun(@times, amat.input, double(~amat.covered))); 
        viscircles([xx,yy],rr,'Color','w','EnhanceVisibility',false,'Linewidth',0.5); 
        viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
        viscircles([xc,yc],rc,'Color','y','EnhanceVisibility',false); 
        title(sprintf('K: covered %d/%d, W: Top-%d disks,\nB: Top-1 disk, Y: previous disk',nnz(amat.covered),H*W,top))
        subplot(223); imshow(amat.axis); title('A-MAT axes')
        subplot(224); imshow(amat.radius,[]); title('A-MAT radii')
        drawnow;
    end
end
amat.reconstruction = reshape(amat.reconstruction,H,W,numChannels);

%% Visualize results
figure(3); clf;
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title('Reconstructed image');
