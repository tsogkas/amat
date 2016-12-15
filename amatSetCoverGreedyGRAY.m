%% Setup global parameters and preprocess image
R = 30; % #scales
B = 32; % #bins
imgRGB  = rgb2gray(im2double(imresize(imread('/home/tsogkas/datasets/AbstractScenes_v1.1/RenderedScenes/Scene0_9.png'),0.5)));
imgLab  = rgb2labNormalized(imgRGB);
[H,W,C] = size(imgRGB);

% Construct filters, calculate perimeneters and disk areas
filters = cell(R,1); for r=1:R, filters{r} = double(disk(r)); end

% Compute encodings f(D_I(x,y,r)) at every point (mean values in this case)
mlab = imageEncoding(imgLab,filters,'average'); 

% Compute sums of I.^2 and I within r-disks for all r.
img = imgRGB; img2 = img.^2; 
sumI = zeros(H,W,R); sumI2 = zeros(H,W,R);
for r=1:R
    sumI(:,:,r)  = conv2(img, filters{r},'same');
    sumI2(:,:,r) = conv2(img2,filters{r},'same');
end

% A(i,j,r) is the #pixels covered by an r-disk centered at (i,j)
A = ones(H,W,R); for r=1:R, A(:,:,r) = conv2(A(:,:,r),filters{r},'same'); end

% We now have to accumulate the squared errors between all r-disks and all
% possible contained disks.
[x,y] = meshgrid(-R:R,-R:R);
dc = bsxfun(@le,x.^2 + y.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
reconstructionError = zeros(H,W,R);
numNewDisksCovered  = zeros(H,W,R);
for r=1:R 
    m = mlab(:,:,1,r);
    N = zeros(H,W);
    s = zeros(H,W);
    s2= zeros(H,W);    
    dcsubset = dc(:,:,end-r+1:end);
    % for a given r-disk, consider all radii of contained disks and
    % accumulate necesasry quantities
    for i=1:size(dcsubset,3) 
%         reconstructionError(:,:,r) = reconstructionError(:,:,r) + conv2(se(:,:,i),D,'same');
        D = dcsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D))); 
        s2= s2 + conv2(sumI2(:,:,i),D,'same');
        s = s  + conv2(sumI(:,:,i), D,'same');
        N = N  + conv2(A(:,:,i),    D,'same');
        numNewDisksCovered(:,:,r) = numNewDisksCovered(:,:,r) + nnz(D);
    end
    reconstructionError(:,:,r) = s2 + N .* m.^2 - 2*m.*s; % add them up    
    % Fix boundary conditions
    reconstructionError([1:r, end-r+1:end],:,r) = inf;
    reconstructionError(:,[1:r, end-r+1:end],r) = inf;
end
reconstructionError = max(0,reconstructionError);
errorBackup = reconstructionError;

% Compute maximality scores using the mean value consensus. Must define
% maximality scores for the disks whose internal part is inside the image
% but the outside ring crosses the image boundary.
A = ones(H,W,R); 
for r=1:R
    dr = ceil(r/10);
    A(:,:,r) = conv2(A(:,:,r),double(ring(r,r+dr)),'same'); 
end
maximalityError = zeros(H,W,R);
mlab2 = mlab.^2;
for r=1:R
%     dr = ceil(r/(2+sqrt(6))); % dA >= 0.5A(r)
    dr = ceil(r/10);
    mask = double(ring(r,r+dr));  % Create masks of outer rings
    maximalityError(:,:,r) = 1-abs(conv2(img,mask,'same') - A(:,:,r).* mlab(:,:,1,r))./A(:,:,r);
%     maximalityError(:,:,r) = 1-(conv2(img2,mask,'same') + ...
%         A(:,:,r).*mlab2(:,:,1,r) - 2 .* mlab(:,:,1,r) .* conv2(img,mask,'same'))./A(:,:,r);
%     maximalityError(:,:,r) = (conv2(imgRGB.^2,mask,'same') + ...
%         nnz(mask)*enc2(:,:,1,r) - 2 .* enc(:,:,1,r) .* conv2(imgRGB,mask,'same'))/nnz(mask);
end 
maximalityError = min(1,max(0,maximalityError));

% Combine the two types of error
A = ones(H,W,R); for r=1:R, A(:,:,r) = conv2(A(:,:,r),filters{r},'same'); end
combinedError = reconstructionError./A + maximalityError;

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
reconstructionError = errorBackup;
numNewPixelsCovered = ones(H,W,R);
for r=1:R
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), filters{r},'same');
end
T = 0.01;
% diskCostEffective = sqrt(reconstructionError ./ numNewPixelsCovered) + bsxfun(@plus,maximalityError,reshape(T./(1:R),1,1,[]));
diskCostEffective = bsxfun(@plus,reconstructionError./numNewPixelsCovered,reshape(T./(1:R).^2,1,1,[]));
% diskCostEffective = bsxfun(@plus,reconstructionError,reshape(T./(1:R),1,1,[]))./numNewPixelsCovered;
% diskCostEffective = reconstructionError ./ numNewPixelsCovered + T*maximalityError;
% diskCostEffective = reconstructionError + maximalityError;
% Sort costs in ascending order and visualize top disks.
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
    % D is the set of points covered by the selected disk
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    % And newPixelsCovered are the NEW points that are covered
    newPixelsCovered = D & ~amat.covered;
    assert(~all(newPixelsCovered(:)))
    
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
    for r=1:R
        xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
        xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
        numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
            numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),filters{r},'same');
        reconstructionError(yymin:yymax,xxmin:xxmax, r) = ...
            reconstructionError(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(priceMap(yymin:yymax,xxmin:xxmax),filters{r},'same');
    end
    reconstructionError = max(0,reconstructionError);

    
    % precompute #new pixels covered by each disk in the subdomain and mark
    % the centers of disks that are fully contained in newPixelsCovered
%     numPixels = zeros(H,W,R);
%     seMap     = zeros(H,W,R);
%     for r=1:R
%         xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
%         xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
%         numPixels(yymin:yymax,xxmin:xxmax,r) = ...
%             conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),filters{r},'same');
%         numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
%             numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - ...
%             numPixels(yymin:yymax,xxmin:xxmax,r);
%         seMap(:,:,r) = (numPixels(:,:,r) == nnz(filters{r})) .* se(:,:,r);
%     end
%     for r=1:R
% %         xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
% %         xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
% %         numPixels = conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),filters{r},'same');
% %         numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
% %             numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - numPixels;
% %         containedDisksInNewArea = numPixels == nnz(filters{r});
%         % se(i,j,r) = squared error of disk of radius r, centered at (i,j)
%         % re(i,j,r) = sum or all se's of contained disks inside the disk of
%         % radius r, centered at (i,j)
%         dcsubset = dc(:,:,end-r+1:end);
%         for i=1:size(dcsubset,3)
%             D = dcsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D)));
%             reconstructionError(yymin:yymax,xxmin:xxmax, r) = ...
%                 reconstructionError(yymin:yymax,xxmin:xxmax, r) - ...
%                 conv2(seMap(yymin:yymax,xxmin:xxmax,i),D,'same');
%         end
%     end
%     assert(~any(reconstructionError(:) < 0))
    
    % Update errors. NOTE: the diskCost for disks that have
    % been completely covered (e.g. the currently selected disk) will be
    % set to inf or nan, because of the division with numNewPixelsCovered
    % which will be zero (0) for those disks. 
%     diskCostEffective = bsxfun(@plus,reconstructionError./numNewPixelsCovered,reshape(T./(1:R).^2,1,1,[]));
    diskCostEffective = reconstructionError ./ numNewPixelsCovered + T*maximalityError;
    diskCostEffective(numNewPixelsCovered == 0) = inf;
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

%% Visualize results
figure(3); clf;
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title('Reconstructed image');
