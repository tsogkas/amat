function mat = refineMAT(mat)

numGroups = numel(unique(mat.branches > 0));
[H,W,~] = size(mat.input);
branches = zeros(H,W,'uint8');
radius   = zeros(H,W);
axis     = zeros(H*W,C);
for i=1:numGroups
    branchOld = mat.branches == i;
    branchNew = bwmorph(bwmorph(branchOld,'dilate',3),'thin',inf);
    radiusOld = branchOld .* double(mat.radius);
    radiusNew = branchNew .* imdilate(radiusOld,strel('disk',2));
    radius(radiusNew > 0) = radiusNew;
    branches(branchNew > 0) = branchNew;
end
[y,x] = find(branches);
r = round(radius(branches));
enc = reshape(permute(mat.enc,[1 2 4 3]), [], C);
idx = sub2ind([H,W,R], y(:),x(:),r(:));
axis(branches>0,:) = enc(idx,:);
mat.axis = axis;
mat.radius = radius;
mat.branches = branches;
% TODO: update depth

