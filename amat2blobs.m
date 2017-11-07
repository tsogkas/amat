function [blobIndIm, blobBoxes, neighbors] = amat2blobs(img, imageSize)

if ischar(img)
    paths = setPaths();
    mat = load(fullfile(paths.amat.precomputed, 'voc2007', ['amat_' img '.mat']));
    mat = mat.mat;
    if nargin < 2
        imageSize = 2*size(mat.depth);
    end
else
    imgResized = imresize(img, 0.5, 'bilinear');
    imgSmoothed = L0Smoothing(imgResized);
    mat = AMAT(imgSmoothed);
    imageSize = [size(img, 1), size(img, 2)];
end
mat.group();
seg = mat.computeSegmentation();
seg = imresize(seg, imageSize, 'nearest');
[blobIndIm, blobBoxes, neighbors] = segAmat2blobsAndNeighbours(seg);

function [blobIndIm, blobBoxes, neighbors] = segAmat2blobsAndNeighbours(seg)

% Fix uncovered parts of the seg
uncovered = seg == 0;
cc = bwconncomp(uncovered);
maxLabel = max(seg(:));
for i=1:cc.NumObjects
    seg(cc.PixelIdxList{i}) = maxLabel + i;
end

% Blob indices
blobIndIm = seg;
segLabels = unique(seg);
numSegments = numel(segLabels);

% Blob boxes and neighbors
blobBoxes = zeros(numSegments, 4);
neighbors = zeros(numSegments);
for i=1:numSegments
    mask = seg == segLabels(i);
    boundaries = bwmorph(mask, 'dilate', 1) & ~mask;
    neighboringSegments = unique(seg(boundaries));
    neighboringSegments(neighboringSegments==0) = [];
    neighbors(i,neighboringSegments) = 1;
    blobBoxes(i,:) = mask2bbox(mask);
end
neighbors = neighbors + eye(numSegments);
assert(issymmetric(neighbors))
