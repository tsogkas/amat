%% Setup global parameters and preprocess image
numScales = 25;
numBins   = 64;
% img0 = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
img0 = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'), [128 128], 'nearest')); 
[H,W,numChannels] = size(img0);
imgLab = rgb2lab(img0);
imgClustered = clusterImageValues(imgLab, 5); % simplify input
% img = imgLab;
img = img0;

%% Construct filters, calculate perimeneters and disk areas
filters = cell(numScales,1); for r=1:numScales, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% Compute individual encodings for color and texture
f = imageEncoding(img,filters);
% f = imageEncoding(binImage(img,numBins),filters,'hist',numBins);

%% Compute decodings g and reconstruction errors at all points and scales
reconstructionError = imageError(img,f,filters,'rse');

%% Greedy approximation of the weighted set cover problem associated with AMAT
% TODO:
% - Manually compare errors and find out why after some point excessively
%   large disks are preferred than smaller disks with zero-error. 
% - Keep track of disk inclusions using sparse matrix and replace for loops
% - Replace nrms error with ssim.
% - Maybe add an extra regularization term that discourages disks with
%   radius that does not agree with the radii of neighboring/enclosed disks

%% Initializations
amat.input          = img;
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
lambda = 0;
reconstructionError = reshape(reconstructionError, H*W,numScales);
reconstructionError = bsxfun(@plus, reconstructionError, lambda./(1:numScales));
costEffectiveness = reconstructionError./ numNewPixelsCovered;
% Sort costs in ascending order and visualize top disks.
[sortedCosts, indSorted] = sort(costEffectiveness(:),'ascend');
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
    [minCost, indMin] = min(costEffectiveness(:));
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
    % LOND VERSION: We have decided to place a disk centered at (xc,yc)
    % with radius r. This disks covers all points (xd,yd) that lie within
    % radius r. For each point (xd,yd) we must find how many disks it is
    % covered by and subtract 1 from the number of NEW pixels covered by
    % these disks (if they were to be selected as part of the solution).
    % Each (xd,yd) is covered by all disks that have centers at distance <=
    % maxRadius from them.
    [yd,xd] = ind2sub([H,W],find(newPixelsCovered)); % find coordinates of all newly covered pixels inside D
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
        end        
    end
    
    % Update cost effectiveness score
    costEffectiveness = reconstructionError./ numNewPixelsCovered;
    assert(all(isinf(costEffectiveness(sub2ind([H,W],yc,xc), 1:rc))))

    
    % Visualize progress
    if 1
        % Sort costs in ascending order to visualize updated top disks.
        [sortedCosts, indSorted] = sort(costEffectiveness(:),'ascend');
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


%% Run the greedy algorithm
reconstructionError = consensusError;
[x,y] = meshgrid(1:W,1:H);
f = mlab;
while ~all(amat.covered(:))
    % Find the most cost-effective set at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    % D is the set of points covered by the selected disk
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    newPixelsCovered = D & ~amat.covered;
    assert(~all(newPixelsCovered(:)))
    
    % Update AMAT
    reconstructedDisk = repmat(reshape(f(yc,xc,:,rc),[1 C]), [nnz(newPixelsCovered),1]);
    amat.price(newPixelsCovered) = sum(( ...
        amat.reconstruction(newPixelsCovered,:) - reconstructedDisk ).^2,2);
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = reconstructedDisk;
    amat.axis(yc,xc,:) = f(yc,xc,:,rc);
    amat.radius(yc,xc) = rc;

    % Find how many of the newPixelsCovered are covered by other disks in
    % the image and subtract the respective counts from those disks.
    % (conv2 for everything, it's the way to go).
    % TODO: HOWEVER, we may have to limit the domain of the convolutions to
    % further speed up, before conv2 gets unnecessarily slow when
    % nnz(newPixelsCovered) is low.
    [yy,xx] = find(newPixelsCovered);
    xmin = min(xx); xmax = max(xx);
    ymin = min(yy); ymax = max(yy);
    newPixelsCovered = double(newPixelsCovered);
    priceMap = amat.price .* newPixelsCovered;
    amat.reconstruction = reshape(amat.reconstruction,H,W,[]); % reshape for convolutions
    for r=1:R
        xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
        xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
        numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
            numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),double(filters{r}),'same');
        reconstructionError(yymin:yymax,xxmin:xxmax, r) = ...
            reconstructionError(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(priceMap(yymin:yymax,xxmin:xxmax),double(filters{r}),'same');
        % TODO: use conv2 with 'valid' parameter to contain and speed up?
        % We must do this "manually" if we want to use imageEncoding() by
        % using the larger window [xxmin,yymin,xxmax,yymax] to compute the
        % convolutions and then crop the [xmin,ymin,xmax,ymax] part.
        tmp = imageEncoding(amat.reconstruction(yymin:yymax,xxmin:xxmax,:),filters(r));
        f(ymin:ymax,xmin:xmax, :, r) = ...
            tmp((ymin-yymin+1):(ymax-yymin+1),(xmin-xxmin+1):(xmax-xxmin+1), :);
    end
    amat.reconstruction = reshape(amat.reconstruction,H*W,[]); % reshape back
    % We do this to correct precision errors introduced by subtracting the
    % pixel prices from the reconstruction error (these would normally be
    % zero but sometimes they take very small negative values).
    reconstructionError = max(0,reconstructionError); 
    
    % Update errors. NOTE: the diskCost for disks that have
    % been completely covered (e.g. the currently selected disk) will be
    % set to inf or nan, because of the division with numNewPixelsCovered
    % which will be zero (0) for those disks. 
    diskCostEffective = min(1,sqrt(reconstructionError ./ numNewPixelsCovered) ...
        + bsxfun(@plus,maximalityError,reshape(T./(1:R),1,1,[])));
    diskCostEffective(isnan(diskCostEffective)) = inf;
    assert(allvec(numNewPixelsCovered(yc,xc, 1:rc)==0))

    
    % Visualize progress
    if 1
        % Sort costs in ascending order to visualize updated top disks.
        [sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
        [yy,xx,rr] = ind2sub([H,W,R], indSorted(1:top));
        subplot(221); imshow(reshape(amat.input, H,W,[])); 
        viscircles([xc,yc],rc, 'Color','k','EnhanceVisibility',false);
        title('Selected disk');
        subplot(222); imshow(bsxfun(@times, reshape(amat.input,H,W,[]), double(~amat.covered))); 
        viscircles([xx,yy],rr,'Color','w','EnhanceVisibility',false,'Linewidth',0.5); 
        viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
        viscircles([xc,yc],rc,'Color','y','EnhanceVisibility',false); 
        title(sprintf('K: covered %d/%d, W: Top-%d disks,\nB: Top-1 disk, Y: previous disk',nnz(amat.covered),H*W,top))
        subplot(223); imshow(amat.axis); title('A-MAT axes')
        subplot(224); imshow(amat.radius,[]); title('A-MAT radii')
        drawnow;
    end
    disp(nnz(~amat.covered))
end
amat.reconstruction = reshape(amat.reconstruction,H,W,C);
amat.input = reshape(amat.input,H,W,C);

