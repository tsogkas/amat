function rec = mat2reconstruction(mat,rad,depth,scales)
% rad = any(mat,3) .* rad: double, HxW array whos nonzero elements equal mat radii

diskf = cell(1,numel(scales)); 
for r=1:numel(diskf) 
    diskf{r} = double(repmat(disk(scales(r)), [1 1 3])); 
end
rind = containers.Map(scales,1:numel(scales));

rec = zeros(size(mat));
rad = round(rad); % make sure radii are integers
[yc,xc] = find(rad);
for p=1:numel(yc)
    x = xc(p); y = yc(p); r = single(rad(y,x)); c = mat(y,x,:);
    rec((y-r):(y+r),(x-r):(x+r),:) = ...
        rec((y-r):(y+r),(x-r):(x+r),:) + bsxfun(@times, diskf{rind(r)}, c);
end
rec = bsxfun(@rdivide, rec, depth);
