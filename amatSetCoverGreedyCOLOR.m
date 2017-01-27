%% Setup global parameters and preprocess image
R = 40; % #scales
B = 32; % #bins
% imgRGB  = im2double(imresize(imread('/home/tsogkas/datasets/AbstractScenes_v1.1/RenderedScenes/Scene0_9.png'),0.5));
imgRGB  = L0Smoothing(im2double(imresize(imread('/home/tsogkas/code/amat/data/BSDS500/images/train/41004.jpg'),0.5)),3e-2);
imgLab  = rgb2labNormalized(imgRGB);
[H,W,C] = size(imgRGB);

% Construct filters, calculate perimeneters and disk areas
filters = cell(R,1); for r=1:R, filters{r} = double(disk(r)); end

% Compute encodings f(D_I(x,y,r)) at every point (mean values in this case)
mlab = imageEncoding(imgLab,filters,'average'); 

% Compute sums of I.^2 and I within r-disks for all r.
img = imgLab; img2 = img.^2; mlab2 = mlab.^2;
sumI2 = zeros(H,W,C,R);
for c=1:C
    for r=1:R
        sumI2(:,:,c,r) = conv2(img2(:,:,c),filters{r}/nnz(filters{r}),'same');
    end
end

% We now have to accumulate the mean squared errors of all contained r-disks
[x,y] = meshgrid(-R:R,-R:R);
dc = bsxfun(@le,x.^2 + y.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
reconstructionCost = zeros(H,W,C,R);
numNewDisksCovered  = zeros(H,W,R);
for r=1:R
    dcsubset = dc(:,:,end-r+1:end);
    % for a given r-disk, consider all radii of contained disks and
    % accumulate necessary quantities
    for i=1:size(dcsubset,3)
        D = dcsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D)));
        for c=1:C
            reconstructionCost(:,:,c,r) = reconstructionCost(:,:,c,r) + ...
                conv2(sumI2(:,:,c,i),  D,'same') + mlab2(:,:,c,r) *nnz(D) - ...
                conv2(mlab(:,:,c,i), D,'same') .* mlab(:,:,c,r) .* 2;
        end
        numNewDisksCovered(:,:,r) = numNewDisksCovered(:,:,r) + nnz(D);
    end
    % Fix boundary conditions
    reconstructionCost([1:r, end-r+1:end],:,:,r) = inf;
    reconstructionCost(:,[1:r, end-r+1:end],:,r) = inf;
end
reconstructionCost = max(0,reconstructionCost); % account for numerical errors (costs should always be positive)
% Combine costs from different channels
wc = [0.5,0.25,0.25]; % weights for luminance and color channels
reconstructionCost = reconstructionCost(:,:,1,:)*wc(1) + ...
                     reconstructionCost(:,:,2,:)*wc(2) + ...
                     reconstructionCost(:,:,3,:)*wc(3);
reconstructionCost = squeeze(reconstructionCost);                 
costBackup = reconstructionCost;

%% Greedy approximation of the weighted set cover problem associated with AMAT
% Initializations
amat.input          = reshape(imgLab, H*W, C);
amat.reconstruction = reshape(amat.input, H*W, C);
amat.axis           = zeros(H,W,C);
amat.radius         = zeros(H,W);
amat.depth          = zeros(H,W); % #disks points(x,y) is covered by
amat.price          = zeros(H,W); % error contributed by each point
amat.se             = zeros(H,W); % squared error at each point
amat.covered        = false(H,W); 
% Flag corners that are not accessible by our filter set
amat.covered([1,end],1)   = true;
amat.covered(end,[1,end]) = true;

% Error balancing and visualization of top (low-cost) disks
% -------------------------------------------------------------------------
% Disk cost: refers to the cost contributed to the total absolute cost of 
% placing the disk in the image.
% -------------------------------------------------------------------------
% Disk cost effective: refers to the "normalized" cost of the disk, divided
% by the number of *new* pixels it covers when it is added in the solution.
% -------------------------------------------------------------------------
% numNewPixelsCovered: number of NEW pixels that would be covered by each
% disk if it were to be added in the solution. This number can be easily 
% computed using convolutions, while taking into account the fact that 
% some disks exceed image boundaries.
% -------------------------------------------------------------------------
reconstructionCost = costBackup;
numNewPixelsCovered = ones(H,W,R);
for r=1:R
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), filters{r},'same');
end
% T = 0.0001;
% diskCostEffective = sqrt(reconstructionError ./ numNewPixelsCovered) + bsxfun(@plus,maximalityError,reshape(T./(1:R),1,1,[]));
% diskCostEffective = bsxfun(@plus,reconstructionCost./numNewPixelsCovered,reshape(T./(1:R).^2,1,1,[]));
% diskCostEffective = bsxfun(@plus,reconstructionError,reshape(T./(1:R),1,1,[]))./numNewPixelsCovered;
% diskCostEffective = reconstructionError ./ numNewPixelsCovered + T*maximalityError;
% diskCostEffective = reconstructionError + maximalityError;
% Sort costs in ascending order and visualize top disks.

% Weights for convex combination of cost types
wm = 1e-7; % maximality coefficient
wr = 1-wm; % reconstruction coefficient
ws = 1e-3;  % scale fixed cost coefficient 
% Define the cost function used to combine the different cost terms
% reconstructionCost = bsxfun(@plus, reconstructionCost , reshape(ws./(1:R),1,1,[])); cf = @() reconstructionCost ./ numNewPixelsCovered; 
cf = @() bsxfun(@plus, reconstructionCost ./ numNewPixelsCovered, reshape(ws./(1:R),1,1,[])); 
% cf = @() bsxfun(@plus, (wr*reconstructionCost + wm*maximalityCost)./ numNewPixelsCovered, reshape(ws./(1:R), 1,1,[])) ;
diskCostEffective = cf();
[sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
top = 1e2;
[yy,xx,rr] = ind2sub([H,W,R], indSorted(1:top));
figure(1); imshow(reshape(amat.input,H,W,[])); 
viscircles([xx,yy],rr,'Color','k','LineWidth',0.5);
viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
title(sprintf('W: Top-%d disks, B: Top-1 disk',top))

%% Run the greedy algorithm
% [xr,yr] = meshgrid(1:R,1:R);
% containedDisks = bsxfun(@le,xr.^2 + yr.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
% coveringDisks  = bsxfun(@le,xr.^2 + yr.^2,reshape((1:R).^2, 1,1,[]));
% [ycov,xcov,rcov] = ind2sub([2*R+1,2*R+1,R],find(dc));
% ycov = ycov-R-1; xcov = xcov-R-1;
[x,y] = meshgrid(1:W,1:H); f = mlab;
while ~all(amat.covered(:))
    % Find the most cost-effective set at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    if isinf(minCost)
        disp('Minimum cost is inf. Stopping execution...')
        break
    end
        
    % D is the set of points covered by the selected disk
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    % And newPixelsCovered are the NEW points that are covered
    newPixelsCovered = D & ~amat.covered;
    if ~any(newPixelsCovered(:))
        disp('All points in the selected disk have been covered.')
        disp('Stopping execution...')
        break
    end
    
    % Update AMAT
    reconstructedDisk = repmat(reshape(f(yc,xc,:,rc),[1 C]), [nnz(newPixelsCovered),1]);
    amat.se(newPixelsCovered) = sum(( ...
        amat.reconstruction(newPixelsCovered,:) - reconstructedDisk ).^2,2);
    amat.price(newPixelsCovered) = minCost / nnz(newPixelsCovered);
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = reconstructedDisk;
    amat.axis(yc,xc,:) = f(yc,xc,:,rc);
    amat.radius(yc,xc) = rc;

    % Find how many of the newPixelsCovered are covered by other disks in
    % the image and subtract the respective counts from those disks.
    [yy,xx] = find(newPixelsCovered);
    xmin = min(xx); xmax = max(xx);
    ymin = min(yy); ymax = max(yy);
    priceMap = amat.price .* newPixelsCovered;
    newPixelsCovered = double(newPixelsCovered);
    costPerPixel = reconstructionCost ./ numNewPixelsCovered;
    for r=1:R
        xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
        xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
        numPixelsSubtracted = conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),filters{r},'same');
        numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
            numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - numPixelsSubtracted;
        reconstructionCost(yymin:yymax,xxmin:xxmax, r) = ...
            reconstructionCost(yymin:yymax,xxmin:xxmax, r) - ...
            numPixelsSubtracted .* costPerPixel(yymin:yymax,xxmin:xxmax, r);
    end
    % Some pixels are assigned NaN values because of the inf-inf
    % subtraction and since max(0,NaN) = 0, we have to reset them
    % explicitly to inf.
    reconstructionCost(isnan(reconstructionCost)) = inf;
%     reconstructionError = max(0,reconstructionError);
    
    % Update errors. NOTE: the diskCost for disks that have
    % been completely covered (e.g. the currently selected disk) will be
    % set to inf or nan, because of the division with numNewPixelsCovered
    % which will be zero (0) for those disks. 
    diskCostEffective = cf();
    diskCostEffective(numNewPixelsCovered == 0) = inf;
    assert(allvec(numNewPixelsCovered(yc,xc, 1:rc)==0))

    
    % Visualize progress
    if 0
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
    fprintf('%d new pixels covered, %d pixels remaining\n',nnz(newPixelsCovered),nnz(~amat.covered))
end

%% Visualize results
amat.reconstruction = reshape(amat.reconstruction,H,W,C);
amat.input = reshape(amat.input,H,W,C);
figure; clf;
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title(['Reconstructed image (ws=' num2str(ws) ')']);
