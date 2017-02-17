function percentage = computeSetCoverage(s)
diskf = cell(1,200);
for r=1:numel(diskf), diskf{r} = disk(r); end

[s(:).nPixelsCovered] = deal(0);
nPixels = 0;
parfor i=1:numel(s)
    [H,W,K] = size(s(i).pts);
    m = false(H,W);
    for k=1:K
        m(:) = false;
        rc = round(s(i).rad(:,:,k));
        [yc,xc] = find(rc);
        for p=1:numel(yc)
            x = xc(p);
            y = yc(p);
            r = single(rc(y,x));
            m((y-r):(y+r),(x-r):(x+r)) = m((y-r):(y+r),(x-r):(x+r)) | diskf{r};
        end
        nPixels = nPixels + nnz(m);
    end
end
percentage = nPixels / sum(cellfun(@numel, {s(:).pts}));

