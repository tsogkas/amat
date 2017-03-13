function mask = mat2mask(rad,scales)
% mat = pts .* rad: double, HxW array whos nonzero elements equal mat radii

diskf = cell(1,numel(scales)); 
for r=1:numel(diskf), diskf{r} = double(disk(scales(r))); end
rind = containers.Map(scales,1:numel(scales));

mask = zeros(size(rad));
rad  = round(rad); % make sure radii are integers
[yc,xc] = find(rad);
for p=1:numel(yc)
    x = xc(p); y = yc(p); r = single(rad(y,x));
    mask((y-r):(y+r),(x-r):(x+r)) = mask((y-r):(y+r),(x-r):(x+r)) + diskf{rind(r)};
end

