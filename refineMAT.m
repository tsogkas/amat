function mat = refineMAT(mat,method,param)

% Default input arguments
if nargin < 3, param  = 3; end
if nargin < 2, method = 'dilation'; end

% Post-processing function
switch method
    case 'dilation'
        SE = strel('disk',param);
        process = @(x) imdilate(x,SE);
    case 'iso-dilation'
        process = @(x) bwdist(x) <= param;
    case 'skeleton'
        process = @(x) bwmorph(x,'skel',inf);
    case 'afmm-skeleton'
        process = @(x) skeleton(x)>=param;
    otherwise
        error(['Method not supported. Supported methods are:\n' ...
                'dilation, iso-dilation, skeleton, afmm-skeleton.'])
end

% The group labels are already sorted and first label is zero (background)
numBranches = max(mat.branches(:)); 
[H,W,C]   = size(mat.input);
branches  = zeros(H,W);
radius    = zeros(H,W);
for i=1:numBranches
    % Old branch points, radii, and respective cover.
    branchOld = mat.branches == i;
    radiusOld = branchOld .* double(mat.radius);
    cover     = mat2mask(radiusOld, mat.scales)>0;
    % Apply post-processing and thinning to selected branch. 
    % Crop the necessary region for more efficiency.
    % Dilation and iso-dilation are applied on the branch points, whereas
    % skeletonization is applied on the cover mask.
    if strcmp(method,'dilation') || strcmp(method,'iso-dilation')
        branchNew = bwmorph(process(branchOld),'thin',inf);
    else
        branchNew = bwmorph(process(cover),'thin',inf);
    end
    % Compute new radii as distance transform on reconstructed cover.
    radiusNew = bwdist(bwperim(cover)) .* double(branchNew);
    % Find closest radii in the subset of the acceptable scale values.
    valid = radiusNew > 0;
    [~,idx] = min(abs(bsxfun(@minus,radiusNew(valid),mat.scales)),[],2);
    radiusNew(valid) = mat.scales(idx);
    % Assign values in global label and radius map.
    branches(valid) = i;
    radius(valid) = radiusNew(valid);
end
assert(all(branches(branches>0) & radius(branches>0)))
assert(all(ismember(radius(radius>0), mat.scales)))

% Make sure there are no gaps among branch labels
newLabels = unique(branches); newLabels(1) = []; % first group is zero
for i=1:numel(newLabels)
    branches(branches == newLabels(i)) = i;
end 

% Find which pixels have been removed and which have been added
oldpts  = any(mat.axis,3);
newpts  = branches > 0;
removed = oldpts & ~newpts;
added   = newpts & ~oldpts;

% Update depth
% NOTE: there is a discrepancy between
% mat2mask(double(newpts).*radius,mat.scales) and 
% mat.depth + depthAdded - depthRemoved. This is probably because when the
% new radii are changed EVEN FOR THE POINTS THAT ARE NOT REMOVED.
% depthAdded   = mat2mask(radius.*double(added),       mat.scales);
% depthRemoved = mat2mask(mat.radius.*double(removed), mat.scales);
depth = mat2mask(radius,mat.scales);

% Update MAT encodings
[y,x] = find(newpts);
r   = radius(newpts);
R   = numel(mat.scales);
enc = reshape(permute(imageEncoding(rgb2labNormalized(mat.input),mat.scales),[1 2 4 3]), [], C);
idx = sub2ind([H,W,R], y(:),x(:),r(:));
newaxis = reshape(rgb2labNormalized(zeros(H,W,C)),H*W,C);
newaxis(newpts,:) = enc(idx,:); % remember that encodings are in LAB!
newaxis = labNormalized2rgb(reshape(newaxis,H,W,C));

% Update reconstruction
reconstruction = mat2reconstruction(reshape(newaxis,H,W,C), radius, depth, mat.scales);

% Update mat fields
mat.radius   = radius;
mat.branches = branches;
mat.axis     = newaxis;
mat.depth    = depth;
mat.reconstruction = reconstruction;



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