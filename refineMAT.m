function mat = refineMAT(mat)

labels = unique(mat.branches);
numGroups = numel(labels)-1; % first label is zero (background)
[H,W,C] = size(mat.input);
branches = zeros(H,W,'uint8');
radius   = zeros(H,W);
axis     = zeros(H*W,C);
for i=1:numGroups
    branchOld = mat.branches == i;
    branchNew = bwmorph(imdilate(branchOld,strel('disk',3)),'thin',inf);
    idx = find(branchNew);
    branches(idx) = branchNew(idx);    
    radiusOld = branchOld .* double(mat.radius);
    cover = mat2mask(radiusOld);
    radiusNew = round(bwdist(bwperim(cover)));
    radiusNew = radiusNew .* double(branchNew); 
    radius(radiusNew>0) = radiusNew(radiusNew>0);
end

assert(allvec(radius(branches>0)>0))
% Update encodings
R = numel(mat.scales);
[y,x] = find(branches);
r = round(radius(branches>0));
enc = reshape(permute(mat.enc,[1 2 4 3]), [], C);
idx = sub2ind([H,W,R], y(:),x(:),r(:));
axis(branches>0,:) = enc(idx,:);
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

% Update mat fields
mat.depth = mat.depth + depthOffset;
mat.axis = axis;
mat.radius = radius;
mat.branches = branches;
