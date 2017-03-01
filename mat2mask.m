function m = mat2mask(mat,diskf)
% mat = pts .* rad: double, HxW array whos nonzero elements equal mat radii

if nargin < 3,
    diskf = cell(1,40); for r=1:numel(diskf), diskf{r} = double(disk(r)); end
end

m = zeros(size(mat));
mat = round(mat);
[yc,xc] = find(mat);
for p=1:numel(yc)
    x = xc(p);
    y = yc(p);
    r = single(mat(y,x));
%     m((y-r):(y+r),(x-r):(x+r)) = m((y-r):(y+r),(x-r):(x+r)) | diskf{r};
    m((y-r):(y+r),(x-r):(x+r)) = m((y-r):(y+r),(x-r):(x+r)) + diskf{r};
end

