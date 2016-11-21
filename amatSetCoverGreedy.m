%% Setup global parameters and preprocess image
numScales = 25;
numBins   = 64;
errorType = 'se';
encodingType = 'average';
colorWeights = [];

% imgRGB = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
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
% Initializations
amat.input          = reshape(imgRGB, H*W, numChannels);
amat.reconstruction = reshape(amat.input, H*W, numChannels);
amat.axis           = zeros(H,W,numChannels);
amat.radius         = zeros(H,W);
amat.depth          = zeros(H,W); % #disks points(x,y) is covered by
amat.covered        = false(H,W);
amat.price          = inf(H,W);   % error contributed by each point


%% Create sparse adjacency matrix for disk inclusions.
% A is a H*W x H*W x R array, with A(i,j,r) == 1 iff the disk that is
% centered at pixel i and has radius r, covers pixel j.
[x,y] = meshgrid(1:W,1:H); 
xx2 = bsxfun(@minus,x(:),x(:)').^2;
yy2 = bsxfun(@minus,y(:),y(:)').^2;
coveredByDisk = false(H*W,H*W,numScales);
for r=1:numScales
    coveredByDisk(:,:,r) = xx2 + yy2 <= r^2;
end
clear xx2 yy2

% Easy way to compute the number of NEW pixels that will be covered by each
% disk if it is added in the solution, taking into account the fact that
% larger disks exceed image boundaries.
numNewPixelsCovered = ones(H,W,numScales);
for r=1:numScales
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), double(filters{r}),'same');
end
numNewPixelsCovered = reshape(numNewPixelsCovered ,H*W,numScales);


%% Error balancing and visualization of top (low-cost) disks
% -------------------------------------------------------------------------
% Disk cost: refers to the cost contributed to the total absolute cost of 
% placing the disk in the image.
% -------------------------------------------------------------------------
% Disk cost effective: refers to the "normalized" cost of the disk, divided
% by the number of *new* pixels it covers when it is added in the solution.
% -------------------------------------------------------------------------
% Before running the algorithm we must add a regularization term, that
% favors larger radii. This is necessary, in order to select "maximal"
% disks. This term must be *added* (not a multiplicative factor), to handle 
% cases when the reconstruction error is 0.
% We consider a threshold T. If the increase in the relative error
% (e_{i+1}/N_{i+1}-e_i/N_i)/(e_i/N_i) < T then the additive factor must be 
% such that it e_{i+1} + f_{i+1} < e_i + f_i. This term must be added to 
% thcost BEFORE normalization, and only in the beginning (no need to
% regularize in the next iterations).
% -------------------------------------------------------------------------
% The final recursive formula for the regularization term is:
% REG(r_R) = 0; REG(r_i) = T*e_i + (N_i / N_{i+1})*REG_{i+1}, for i=1,...,R-1
T = 0.0001;
diskCost = reshape(reconstructionError, H*W,numScales);
freshaped = reshape(f, [],numChannels,numScales);
REG = [squeeze(mean(diff(freshaped,1,3).^2 ,2)), zeros(H*W,1)]; 
REG = bsxfun(@plus, REG, (numScales:-1:1)*T);
diskCostEffective = diskCost ./ numNewPixelsCovered + REG;
% Sort costs in ascending order and visualize top disks.
[sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
top = 1e2;
[yy,xx,rr] = ind2sub([H,W,numScales], indSorted(1:top));
figure(1); imshow(amat.input); 
viscircles([xx,yy],rr,'Color','k','LineWidth',0.5);
viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
title(sprintf('W: Top-%d disks, B: Top-1 disk',top))


%% Run the greedy algorithm
[x,y] = meshgrid(1:W,1:H);
while ~all(amat.covered(:))
    % Find the most cost-effective set at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    % D is the set of points covered by the selected disk
    [yc,xc,rc] = ind2sub([H,W,numScales], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    newPixelsCovered = D & ~amat.covered;
    
    % Update AMAT ()
    reconstructedDisk = repmat(reshape(f(yc,xc,:,rc),[1 3]), [nnz(newPixelsCovered),1]);
    amat.price(newPixelsCovered) = mean(( ...
        amat.reconstruction(newPixelsCovered,:) - reconstructedDisk ).^2,2);
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = reconstructedDisk;
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
    [yd,xd] = ind2sub([H,W],find(newPixelsCovered)); % #pixels x 1
    xx2 = bsxfun(@minus,xd,x(:)').^2; % #pixels x H*W
    yy2 = bsxfun(@minus,yd,y(:)').^2; % #pixels x H*W
    % isCoveredByDisk: #pixels x H*W x numScales
    % isCoveredByDisk(i,j,r): pixel i is covered by a disk centered at
    % pixel j, with radius r. 
    % sum(isCoveredByDisk,1): #pixels covered by pixel j at radius r
    isCoveredByDisk = bsxfun(@le,xx2+yy2,reshape((1:numScales).^2,1,1,[])); 
    numNewPixelsCovered = numNewPixelsCovered - squeeze(sum(isCoveredByDisk));
    % Remove error caused by pixels that are covered by the selected disk.
    
    
    
    
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
            diskCost(inds,r) = diskCost(inds,r) - ...
                mean(bsxfun(@minus, freshaped(inds,:,r), ... % subtract old cost for pixel
                amat.input(sub2ind([H,W],yd(i),xd(i)),:)).^2, 2);
%             freshaped(inds,:,r) = freshaped(inds,:,r) - 
        end        
    end
    
    % Update cost effectiveness score
    diskCostEffective = diskCost./ numNewPixelsCovered;
    assert(allvec(numNewPixelsCovered(sub2ind([H,W],yc,xc), 1:rc)==0))

    
    % Visualize progress
    if 1
        % Sort costs in ascending order to visualize updated top disks.
        [sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
        [yy,xx,rr] = ind2sub([H,W,numScales], indSorted(1:top));
        subplot(221); imshow(reshape(amat.input, H,W)); 
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
