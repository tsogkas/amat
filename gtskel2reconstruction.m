function rec = gtskel2reconstruction(img,pts,rad)
scales = 2:41;
diskf = cell(1,numel(scales)); 
for r=1:numel(diskf), diskf{r} = double(disk(scales(r))); end
rind = containers.Map(scales,1:numel(scales));

[H,W,C]   = size(img);
img       = im2double(img);
numSkels  = size(pts,3);
nnzPixels = zeros(numSkels,1);
rec       = zeros(H,W,C,numSkels);
for i=1:numSkels
    skel = bwthin(pts(:,:,i));
    skel(border(skel,scales(1))) = 0;
    rads = double(skel) .* double(rad(:,:,i));
    % Adjust scales and skeletons
    valid = rads > 0;
    [~,idxScales] = min(abs(bsxfun(@minus,rads(valid),scales)),[],2);
    rads(valid) = scales(idxScales);
    % Compute reconstruction
    depth = zeros(size(skel));
    [yc,xc] = find(rads);
    for p=1:numel(yc)
        x = xc(p); y = yc(p); r = rads(y,x); 
        % Because of resizing, r sometimes gets out of bound
        r = min(r,min(x-1,y-1)); r = min(r, min(H-y, W-x));
        dfilt = diskf{rind(r)};
        patch = bsxfun(@times, img(y,x,:), dfilt);
        mval  = sum(sum(patch))/nnz(dfilt);
        patch = bsxfun(@times, mval, dfilt);
        rec((y-r):(y+r),(x-r):(x+r),:,i) = ...
            rec((y-r):(y+r),(x-r):(x+r),:,i) + patch;
        depth((y-r):(y+r),(x-r):(x+r)) = ...
            depth((y-r):(y+r),(x-r):(x+r)) + dfilt;
    end
    rec(:,:,:,i) = bsxfun(@rdivide, rec(:,:,:,i), depth); 
    rec(:,:,:,i) = reshape(inpaint_nans(rec(:,:,:,i)), H,W,[]);
    nnzPixels(i) = nnz(skel);
end