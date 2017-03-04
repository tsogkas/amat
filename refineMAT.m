function mat = refineMAT(mat)

% The group labels are already sorted and first label is zero (background)
numBranches = max(mat.branches(:)); 
[H,W,C]   = size(mat.input);
branches  = zeros(H,W);
radius    = zeros(H,W);
mataxes   = reshape(rgb2labNormalized(zeros(H,W,C)),H*W,C);
SE        = ones(2);
for i=1:numBranches
    branchOld = mat.branches == i;
    branchNew = bwmorph(imdilate(branchOld,SE),'thin',inf);
    radiusOld = branchOld .* double(mat.radius);
    cover     = mat2mask(radiusOld);
    radiusNew = round(bwdist(bwperim(cover)));
    radiusNew = radiusNew .* double(branchNew);
    branchNew(branchNew > 0 & radiusNew == 0) = 0;
    branches(branchNew) = i;
    radius(radiusNew>0) = radiusNew(radiusNew>0);
end

% Re-adjust labels
newGroups = unique(branches); newGroups(1) = []; % first group is zero
for i=1:numel(newGroups)
    branches(branches == newGroups(i)) = i;
end

assert(all(radius(branches>0)>0))
% Update encodings
R = numel(mat.scales);
[y,x] = find(branches);
r = round(radius(branches>0));
enc = reshape(permute(mat.enc,[1 2 4 3]), [], C);
idx = sub2ind([H,W,R], y(:),x(:),r(:));
mataxes(branches>0,:) = enc(idx,:);

% Update depth
removed = any(mat.axis,3) & ~(branches > 0);
added   = (branches > 0) & ~any(mat.axis,3);
depthOffset = zeros(H,W);
[xx,yy] = meshgrid(1:W,1:H);

% First subtract depth of removed pixels
[y,x] = find(removed);
for i=1:numel(y)
    r = mat.radius(y(i),x(i)); % radii of old disks
    xmin = max(1,x(i)-r); xmax = min(W,x(i)+r);
    ymin = max(1,y(i)-r); ymax = min(H,y(i)+r);
    depthOffset(ymin:ymax,xmin:xmax) = depthOffset(ymin:ymax,xmin:xmax) - ...
        double(((xx(ymin:ymax,xmin:xmax)-x(i)).^2 + ...
         (yy(ymin:ymax,xmin:xmax)-y(i)).^2) <= r^2);
end

% Then add depth of added pixels
[y,x] = find(added);
for i=1:numel(y)
    r = radius(y(i),x(i)); % radii of new disks
    xmin = max(1,x(i)-r); xmax = min(W,x(i)+r);
    ymin = max(1,y(i)-r); ymax = min(H,y(i)+r);
    depthOffset(ymin:ymax,xmin:xmax) = depthOffset(ymin:ymax,xmin:xmax) + ...
        double(((xx(ymin:ymax,xmin:xmax)-x(i)).^2 + ...
         (yy(ymin:ymax,xmin:xmax)-y(i)).^2) <= r^2);
end

% Reconstruct the reconstruction
reconstruction = zeros(H*W,C);


% Update mat fields
mat.depth = mat.depth + depthOffset;
mat.axis = labNormalized2rgb(reshape(mataxes,H,W,C));
mat.radius = radius;
mat.branches = branches;
