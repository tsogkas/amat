function mat = amat(img,scales,ws,vistop)

if nargin < 2, scales = 40; end
if nargin < 3, ws = 1e-3; end
if nargin < 4, vistop = 0; end
if isscalar(scales)
    R = scales;
    scales = 1:R;
elseif isvector(scales)
    R = numel(scales);
else
    error('''scales'' must be a vector of disk radii or a scalar (#scales)')
end

% Convert to CIE La*b* color space
img = rgb2labNormalized(im2double(img));

% Compute encodings f(D_I(x,y,r)) at every point (mean values in this case)
enc = imageEncoding(img,scales);

% Compute disk cost based on image reconstruction heuristic
% profile on;
diskCost = reconstructionCost(enc,scales);

% Compute reedy approximation of the weighted set cover for the AMAT.
mat = setCover(img,enc,diskCost,scales,ws,vistop);
% profile off; profile viewer;
end

% -------------------------------------------------------------------------
function mat = setCover(img,enc,diskCost,scales,ws,vistop)
% -------------------------------------------------------------------------
% Greedy approximation of the weighted set cover problem.
% -------------------------------------------------------------------------
% - Disk cost: cost incurred by selecting ena r-disk, centered at (i,j).
% - numNewPixelsCovered: number of NEW pixels covered by a selected disk.
% - Cost per pixel: diskCost / numNewPixelsCovered.
% - Disk cost effective: adjusted normalized cost: costPerPixel + scaleTerm
%       where scaleTerm is a term that favors selecting disks of larger
%       radii. Such a term is necessary, to resolve selection of disks in
%       the case where diskCost is zero for more than on radii.
% 
% TODO: is there a way to first sort scores and then pick the next one in
%       queue, to avoid min(diskCostEffective(:)) in each iteration?
% TODO: check why numNewPixelsCovered and nnz(newPixelsCovered) are not
%       always equal.

% Initializations
[H,W,C,R]          = size(enc);
mat.input          = reshape(img, H*W, C);
mat.reconstruction = reshape(mat.input, H*W, C);
mat.axis           = zeros(H,W,C);
mat.radius         = zeros(H,W);
mat.depth          = zeros(H,W); % #disks points(x,y) is covered by
mat.price          = zeros(H,W); % error contributed by each point
mat.se             = zeros(H,W); % squared error at each point
mat.covered        = false(H,W); 
% Flag corners that are not accessible by our filter set
mat.covered([1,end],1)   = true;
mat.covered(end,[1,end]) = true;
mat.ws = ws; % weight for scale-dependent cost term.

% Compute how many pixels are covered be each r-disk.
filters = cell(1,R); for r=1:R, filters{r} = double(disk(scales(r))); end
numNewPixelsCovered = ones(H,W,R);
for r=1:R
    numNewPixelsCovered(:,:,r) = ...
        conv2(numNewPixelsCovered(:,:,r), filters{r},'same');
end

% Add scale-dependent cost term to favor the selection of larger disks.
costPerPixel = diskCost ./ numNewPixelsCovered;
diskCostEffective = bsxfun(@plus, costPerPixel, reshape(ws./scales,1,1,[]));

[x,y] = meshgrid(1:W,1:H);
% while ~all(mat.covered(:))
for i=1:1000
    % Find the most cost-effective disk at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    if isinf(minCost), break; end
        
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= scales(rc)^2;   % points covered by the selected disk
    newPixelsCovered  = D & ~mat.covered;        % NEW pixels that are covered by D
    
    % Update MAT 
    mat = updateMAT(mat);
    if ~any(newPixelsCovered(:)), break; end
    
    % Update costs 
    [yy,xx] = find(newPixelsCovered);
    xmin = min(xx); xmax = max(xx);
    ymin = min(yy); ymax = max(yy);
    newPixelsCovered = double(newPixelsCovered);
    for r=1:R
        scale = scales(r);
        x1 = max(xmin-scale,1); y1 = max(ymin-scale,1);
        x2 = min(xmax+scale,W); y2 = min(ymax+scale,H);
        % Find how many of the newPixelsCovered are covered by other disks.
        numPixelsSubtracted = ...
            conv2(newPixelsCovered(y1:y2,x1:x2),filters{r},'same');
        % and subtract the respective counts from those disks.
        numNewPixelsCovered(y1:y2,x1:x2, r) = ...
            numNewPixelsCovered(y1:y2,x1:x2, r) - numPixelsSubtracted;
        % update diskCost, costPerPixel, and diskCostEfficiency *only* for
        % the locations that have been affected, for efficiency.
        diskCost(y1:y2,x1:x2, r) = max(0, diskCost(y1:y2,x1:x2, r) - ...
            numPixelsSubtracted .* costPerPixel(y1:y2,x1:x2, r));
        costPerPixel(y1:y2,x1:x2, r) = diskCost(y1:y2,x1:x2, r) ./ ...
            max(eps,numNewPixelsCovered(y1:y2,x1:x2, r)) + ... % avoid 0/0
            inf*(numNewPixelsCovered(y1:y2,x1:x2, r) == 0);    % x/0 = inf
        diskCostEffective(y1:y2,x1:x2, r) = ...
            costPerPixel(y1:y2,x1:x2, r) + ws/scales(r);
    end
    assert(allvec(numNewPixelsCovered(yc,xc, 1:rc)==0))
    
    if vistop, visualizeProgress(mat,diskCostEffective); end
    fprintf('%d new pixels covered, %d pixels remaining\n',...
        nnz(newPixelsCovered),nnz(~mat.covered))
end
mat.reconstruction = reshape(mat.reconstruction,H,W,C);
mat.input = reshape(mat.input,H,W,C);
mat.visualize = @()visualize(mat);

% -------------------------------------------------------------------------
function mat = updateMAT(mat)
% -------------------------------------------------------------------------
reconstructedDisk = ...
    repmat(reshape(enc(yc,xc,:,rc),[1 C]), [nnz(newPixelsCovered),1]);
mat.se(newPixelsCovered) = sum(( ...
    mat.reconstruction(newPixelsCovered,:) - reconstructedDisk ).^2,2);
mat.price(newPixelsCovered) = minCost / nnz(newPixelsCovered);
mat.covered(D) = true;
mat.depth(D) = mat.depth(D) + 1;
mat.reconstruction(newPixelsCovered,:) = reconstructedDisk;
mat.axis(yc,xc,:) = enc(yc,xc,:,rc);
mat.radius(yc,xc) = scales(rc);
end

% -------------------------------------------------------------------------
function visualizeProgress(mat,diskCost)
% -------------------------------------------------------------------------
% Sort costs in ascending order to visualize updated top disks.
[~, indSorted] = sort(diskCost(:),'ascend');
[yy,xx,rr] = ind2sub([H,W,R], indSorted(1:vistop));
subplot(221); imshow(reshape(mat.input, H,W,[]));
viscircles([xc,yc],rc, 'Color','k','EnhanceVisibility',false); title('Selected disk');
subplot(222); imshow(bsxfun(@times, reshape(mat.input,H,W,[]), double(~mat.covered)));
viscircles([xx,yy],rr,'Color','w','EnhanceVisibility',false,'Linewidth',0.5);
viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false);
viscircles([xc,yc],rc,'Color','y','EnhanceVisibility',false);
title(sprintf('K: covered %d/%d, W: Top-%d disks,\nB: Top-1 disk, Y: previous disk',nnz(mat.covered),H*W,vistop))
subplot(223); imshow(mat.axis); title('A-MAT axes')
subplot(224); imshow(mat.radius,[]); title('A-MAT radii')
drawnow;
end

% -------------------------------------------------------------------------
function visualize(mat)
% -------------------------------------------------------------------------
figure; clf;
subplot(221); imshow(mat.axis);             title('Medial axes');
subplot(222); imshow(mat.radius,[]);        title('Radii');
subplot(223); imshow(mat.input);            title('Original image');
subplot(224); imshow(mat.reconstruction);   title('Reconstructed image');
end

end


% -------------------------------------------------------------------------
function cost = reconstructionCost(enc,scales)
% -------------------------------------------------------------------------
% img: HxWxC input image (in the Lab color space). 
% enc: HxWxCxR appearance encodings for all r-disks (mean RGB values by default).
% filters: disk-shaped filters.
% 
% This function computes a heuristic that represents the ability to
% reconstruct a disk-shaped part of the input image, using the mean RGB
% values computed over the same area. Intuitively, the idea behind this
% heuristic is the following: 
% In order to accurately reconstruct an image disk of radius r, using its
% mean RGB values, we must also be able to reconstruct *every* fully
% contained disk of radius r' < r (uniformity criterion).
% 
% Terminology: an r-disk is a disk of radius = r.

% Create grid of relative positions of disk centers for all scales. 
[H,W,C,R] = size(enc);
Rmax = scales(end); [x,y] = meshgrid(-Rmax:Rmax,-Rmax:Rmax);

% cd is a HxWxR logical array, composed of 2D binary masks. 
% cd(i,j,k) is true iff a "scales(end-1+k)-disk", centered at (i,j)
% is *fully contained* within a "scales(end-1+n)-disk", centered at (0,0),
% where n > k.
% In other words, each nonzero point in cd(:,:,k+1:end) corresponds to a
% different disk, contained in the disk represented by the binary mask in
% cd(:,:,k).
% TODO: Should I add the scales from scales(1):-1:0 too?
cd = bsxfun(@le, x.^2+y.^2, reshape([scales(end-1:-1:1) scales(1)-1].^2, 1,1,[]));
enc2  = enc.^2;
encx2 = 2*enc;
cost = zeros(H,W,C,R);
for c=1:C
    for r=1:R
        cdsubset = cd(:,:,end-r+1:end);
        % for a given r-disk, consider all contained disks and accumulate
        % M_R = sum(I_R)/D_R; M_ri = sum(I_ri)/R_ri;
        % Cost = sum((M_R-M_ri)^2) for all enclosed ri-disks.
        % Cost = sum( M_R^2 + M_ri^2 - 2*M_R*M_ri ) = ...
        % D_R*enc2 + conv2(enc2) + encx2 .* conv2(enc)
        encx2t = encx2(:,:,c,r);
        for i=1:size(cdsubset,3)
            D = cdsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D)));
            cost(:,:,c,r) = cost(:,:,c,r) + ...
                conv2(enc2(:,:,c,i), D,'same') - ...
                conv2(enc(:,:,c,i),  D,'same').* encx2t;
        end
        cost(:,:,c,r) = cost(:,:,c,r) + enc2(:,:,c,r) * nnz(cdsubset);
    end
end
tmp = cost;

filters = cell(1,R); 
for r=2:R, filters{r} = double(circle(scales(r-1))); end
filters{1} = disk(scales(1)-1);
enccsum = cumsum(enc,4);
enc2csum= cumsum(enc2,4);
for c=1:C
    for r=1:R
    end
end

% Accumulate the mean squared errors of all contained r-disks.
cost = zeros(H,W,C,R);
for c=1:C
    for r=1:R
        cdsubset = cd(:,:,end-r+1:end);
        % for a given r-disk, consider all contained disks and accumulate
        % M_R = sum(I_R)/D_R; M_ri = sum(I_ri)/R_ri;
        % Cost = sum((M_R-M_ri)^2) for all enclosed ri-disks.
        sumMri = zeros(H,W);
        for i=1:size(cdsubset,3)
            D = cdsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D)));
            sumMri = sumMri + conv2(enc(:,:,c,i),  D,'same');
        end
        cost(:,:,c,r) = (enc(:,:,c,r)*nnz(cdsubset) - sumMri).^2;
    end
end


% Fix boundary conditions. Setting r-borders to a very big cost helps us 
% avoid selecting disks that cross the image boundaries.
% We do not use Inf to avoid complications in the greedy set cover 
% algorithm, caused by inf-inf subtractions and inf/inf divisions.
% Also, keep in mind that max(0,NaN) = 0.
BIG = 1e30;
for r=1:R    
    cost([1:r, end-r+1:end],:,:,r) = BIG;
    cost(:,[1:r, end-r+1:end],:,r) = BIG;
end

% Sometimes due to numerical errors, cost are slightly negative
cost = max(0,cost); 

% Combine costs from different channels
if C > 1
    wc = [0.5,0.25,0.25]; % weights for luminance and color channels
    cost = cost(:,:,1,:)*wc(1) + cost(:,:,2,:)*wc(2) + cost(:,:,3,:)*wc(3);
    cost = squeeze(cost);                 
end
end

% -------------------------------------------------------------------------
function enc = imageEncoding(img,scales)
% -------------------------------------------------------------------------
% Fast version of imageEncoding, using convolutions with circles + cumsum
% instead of convolutions with disks. 
[H,W,C] = size(img); R = numel(scales);
filters = cell(1,R); for r=1:R, filters{r} = double(circle(scales(r))); end
filters{1}(2,2) = 1; % if we use circle filters, center pixel is left out
enc = zeros(H,W,C,R);
for c=1:C
    for r=1:R
        enc(:,:,c,r) = conv2(img(:,:,c),filters{r},'same');
    end
end
enc = cumsum(enc,4);
areas = cumsum(cellfun(@nnz,filters));
enc = bsxfun(@rdivide, enc, reshape(areas,1,1,1,[]));
end