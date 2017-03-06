function mat = amat(img,scales,ws,vistop)

if nargin < 2, scales = 2:41; end
if nargin < 3, ws = 1e-4; end
if nargin < 4, vistop = 0; end
if isscalar(scales)
    scales = 2:(scales+1);
elseif ~isvector(scales)
    error('''scales'' must be a vector of disk radii or a scalar (#scales)')
end

% Convert to CIE La*b* color space
img = rgb2labNormalized(im2double(img));

% Compute encodings f(D_I(x,y,r)) at every point (mean values in this case)
enc = imageEncoding(img,scales);

% Compute disk cost based on image reconstruction heuristic
diskCost = reconstructionCost(enc,scales);

% Compute greedy approximation of the weighted set cover for the AMAT.
% profile on;
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

% Initializations
[H,W,C,R]          = size(enc);
zeroLabNormalized  = rgb2labNormalized(zeros(H,W,C,'single'));
mat.input          = reshape(img, H*W, C);
mat.reconstruction = reshape(zeroLabNormalized,H*W,C);
mat.axis           = zeroLabNormalized;
mat.radius         = zeros(H,W,'single');
mat.depth          = zeros(H,W,'single'); % #disks points(x,y) is covered by
mat.price          = zeros(H,W,'single'); % error contributed by each point
mat.covered        = false(H,W); 
% Flag border pixels that cannot be accessed by filters.
r = scales(1);
mat.covered([1:r,end-r+1:end], [1,end]) = true;
mat.covered([1,end], [1:r,end-r+1:end]) = true;
mat.ws = ws; % weight for scale-dependent cost term.
mat.scales = scales;
BIG = 1e30;

% Compute how many pixels are covered be each r-disk.
filters = cell(1,R); for r=1:R, filters{r} = double(disk(scales(r))); end
diskAreas = cellfun(@nnz,filters);
numNewPixelsCovered = repmat(reshape(diskAreas,1,1,[]), [H,W]);

% Add scale-dependent cost term to favor the selection of larger disks.
costPerPixel = diskCost ./ numNewPixelsCovered;
diskCostEffective = bsxfun(@plus, costPerPixel, reshape(ws./scales,1,1,[]));

% Print remaining pixels to be covered in these points
printBreakPoints = floor((4:-1:1).*(H*W/5));

fprintf('Pixels remaining: ');
[x,y] = meshgrid(1:W,1:H);
while ~all(mat.covered(:))
    % Find the most cost-effective disk at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    if isinf(minCost), 
        warning('Stopping: selected disk has infinite cost.')
        break; 
    end
        
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    D = (x-xc).^2 + (y-yc).^2 <= scales(rc)^2; % points covered by the selected disk
    newPixelsCovered  = D & ~mat.covered;      % NEW pixels that are covered by D
    if ~any(newPixelsCovered(:))
        warning('Stopping: selected disk covers zero (0) new pixels.')
        break; 
    end
    
    % Update MAT 
    mat = updateMAT(mat);
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
        diskCost(y1:y2,x1:x2, r) = diskCost(y1:y2,x1:x2, r) - ...
            numPixelsSubtracted .* costPerPixel(y1:y2,x1:x2, r);
        costPerPixel(y1:y2,x1:x2, r) = diskCost(y1:y2,x1:x2, r) ./ ...
            max(eps,numNewPixelsCovered(y1:y2,x1:x2, r)) + ... % avoid 0/0
            BIG*(numNewPixelsCovered(y1:y2,x1:x2, r) == 0);    % x/0 = inf
        diskCostEffective(y1:y2,x1:x2, r) = ...
            costPerPixel(y1:y2,x1:x2, r) + ws/scales(r);
    end
    % Make sure the same point is not selected again
    diskCost(yc,xc,:) = BIG; diskCostEffective(yc,xc,:) = BIG;
    
    
    if vistop, visualizeProgress(mat,diskCostEffective); end
    if ~isempty(printBreakPoints) && nnz(~mat.covered) < printBreakPoints(1)
        fprintf('%d...',printBreakPoints(1))
        printBreakPoints(1) = [];
    end
end
fprintf('\n')
mat.input = labNormalized2rgb(reshape(mat.input,H,W,C));
mat.axis  = labNormalized2rgb(mat.axis);
mat.reconstruction = mat2reconstruction(mat.axis,mat.radius,mat.depth,mat.scales);
mat.visualize = @()visualize(mat);

% -------------------------------------------------------------------------
function mat = updateMAT(mat)
% -------------------------------------------------------------------------
mat.price(newPixelsCovered) = minCost / numNewPixelsCovered(yc,xc,rc);
mat.covered(newPixelsCovered) = true;
mat.depth(D) = mat.depth(D) + 1;
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
subplot(221); imshow(mat.axis);             title('Medial axes');
subplot(222); imshow(mat.radius,[]);        title('Radii');
subplot(223); imshow(mat.input);            title('Original image');
subplot(224); imshow(mat.reconstruction);   title('Reconstructed image');
end

end


% -------------------------------------------------------------------------
function enc = imageEncoding(img,scales)
% -------------------------------------------------------------------------
% Fast version of imageEncoding, using convolutions with circles + cumsum
% instead of convolutions with disks. 
[H,W,C] = size(img); R = numel(scales);
filters = cell(1,R); filters{1} = double(disk(scales(1))); 
for r=2:R, filters{r} = double(circle(scales(r))); end
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
% 
% The heuristic we use sums all square errors between the encoding of an 
% r-disk centered at a point (i,j) and the encodings of all FULLY CONTAINED
% disks. Written in a simplified mathematical form, for a given r_k-disk:
% M_rk = sum_i(I_rk)/D_rk; M_ri = sum_i(I_ri)/R_ri;
% Cost = sum((M_rk-M_ri)^2) for all enclosed ri-disks.
% Cost = sum( M_rk^2 + M_ri^2 - 2*M_rk*M_ri ) = ...
% D_rk*enc2 + conv2(enc2) + 2 .* enc .* conv2(enc)
% Given an r-disk, filters(r-i+1) is a mask that marks the centers of all
% contained i-disks.

[H,W,C,R]  = size(enc);
filters    = cell(1,R); 
filters{1} = double(disk(scales(1)-1));
for r=2:R, filters{r} = double(circle(scales(r-1))); end
% Precompute necessary quantitities. We use circular filters applied on
% cumulative sums instead of disk filters, for efficiency. 
enc2     = enc.^2;
enccsum  = cumsum(enc,4);
enc2csum = cumsum(enc2,4);
nnzcd    = cumsum(cumsum(cellfun(@nnz, filters)));

cost = zeros(H,W,C,R);
for c=1:C
    for r=1:R
        sumMri  = zeros(H,W);
        sumMri2 = zeros(H,W);
        for i=1:r
            sumMri  = sumMri  + conv2(enccsum(:,:,c,i), filters{r-i+1},'same');
            sumMri2 = sumMri2 + conv2(enc2csum(:,:,c,i),filters{r-i+1},'same');
        end
        cost(:,:,c,r) = enc2(:,:,c,r)*nnzcd(r) + sumMri2 - 2*enc(:,:,c,r).*sumMri;
    end
end

% Fix boundary conditions. Setting scale(r)-borders to a very big cost 
% helps us avoid selecting disks that cross the image boundaries.
% We do not use Inf to avoid complications in the greedy set cover 
% algorithm, caused by inf-inf subtractions and inf/inf divisions.
% Also, keep in mind that max(0,NaN) = 0.
BIG = 1e30;
for r=1:R
    scale = scales(r);
    cost([1:scale, end-scale+1:end],:,:,r) = BIG;
    cost(:,[1:scale, end-scale+1:end],:,r) = BIG;
end

% Sometimes due to numerical errors, cost are slightly negative. Fix this.
cost = max(0,cost); 

% Combine costs from different channels
if C > 1
    wc = [0.5,0.25,0.25]; % weights for luminance and color channels
    cost = cost(:,:,1,:)*wc(1) + cost(:,:,2,:)*wc(2) + cost(:,:,3,:)*wc(3);
    cost = squeeze(cost);                 
end
end
