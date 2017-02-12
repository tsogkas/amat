function m = mat2mask(pts,rad,diskf)

if nargin < 3,
    diskf = cell(1,200); for r=1:numel(diskf), diskf{r} = disk(r); end
end

[H,W] = size(pts);
m = false(H,W);
[yc,xc] = find(rad);
for pts=1:numel(yc)
    x = single(xc(pts));
    y = single(yc(pts));
    r = single(rad(y,x));
    d = m((y-r):(y+r),(x-r):(x+r));
    d(diskf{r}) = true;
    m((y-r):(y+r),(x-r):(x+r)) = d;
end

