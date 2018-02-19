function boxes  = amatBoxes(img, imgId)

% If an image id string is not provided, we will have to compute the AMAT
% from scratch
if nargin < 2 || isempty(imgId)
    imgId = img;
end

% Here you specify which similarity functions to use in merging
simFunctionHandles = {@SSSimColourTextureSizeFillOrig, ...
    @SSSimTextureSizeFill, ...
    @SSSimBoxFillOrig, ...
    @SSSimSize};

seg = amat2seg(imgId, [size(img,1), size(img,2)]);
[blobIndIm, blobBoxes, neighbors] = segAmat2blobsAndNeighbours(seg);
[boxes, ~, priority] = seg2hierarchicalGrouping( im2double(img), ...
    blobIndIm, blobBoxes, neighbors, simFunctionHandles);

% Do pseudo random sorting as in the Selective Search paper
priority = priority .* rand(size(priority));
[~, sortIds] = sort(priority, 'ascend');
boxes = boxes(sortIds,:);

% Filter small and duplicate boxes
minBoxWidth = 20;
boxes = FilterBoxesWidth(boxes, minBoxWidth);
boxes = BoxRemoveDuplicates(boxes);

% -------------------------------------------------------------------------
function seg = amat2seg(img, imageSize)
% -------------------------------------------------------------------------
if ischar(img)
    % Passing the image id string assumes we have precomputed its AMAT.
    paths = setPaths();
    mat = load(fullfile(paths.amat.precomputed, 'voc2007', ['amat_' img '.mat']));
    mat = mat.mat;
    % Resize to a specified dimension, otherwise we assume that the AMAT
    % was computed on an input that has been downsampled by a factor of 2.
    if nargin < 2
        imageSize = 2*size(mat.depth);
    end
else
    % If we are passing the image itself, compute the AMAT from scratch
    % (downsample by a factor of 2 first for efficiency)
    imgResized = imresize(img, 0.5, 'bilinear');
    imgSmoothed = L0Smoothing(imgResized);
    mat = AMAT(imgSmoothed);
    imageSize = [size(img, 1), size(img, 2)];
end
mat.group();
seg = mat.computeSegmentation();
seg = imresize(seg, imageSize, 'nearest');

% -------------------------------------------------------------------------
function [blobIndIm, blobBoxes, neighbors] = segAmat2blobsAndNeighbours(seg)
% -------------------------------------------------------------------------
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
% We use the format used in the Selective Search functions for consistency
blobBoxes = blobBoxes(:,[2,1,4,3]);
neighbors = neighbors + eye(numSegments);
assert(issymmetric(neighbors))

% -------------------------------------------------------------------------
function [boxes, hierarchy, priority] = seg2hierarchicalGrouping(...
    colourIm, blobIndIm, blobBoxes, neighbours, functionHandles)
% -------------------------------------------------------------------------
% Copied part of the Image2HierarchicalGrouping() function, found in the
% original Selective Search code.

numBlobs = size(blobBoxes,1);

% Skip hierarchical grouping if segmentation results in single region only
if numBlobs == 1
    warning('Oversegmentation results in a single region only');
    boxes = blobBoxes;
    hierarchy = [];
    priority = 1; % priority is legacy
    return;
end

%%% Calculate histograms and sizes as prerequisite for grouping procedure

% Get colour histogram
[colourHist, blobSizes] = BlobStructColourHist(blobIndIm, colourIm);

% Get texture histogram
textureHist = BlobStructTextureHist(blobIndIm, colourIm);
% textureHist = BlobStructTextureHistLBP(blobIndIm, colourIm);

% Allocate memory for complete hierarchy.
blobStruct.colourHist = zeros(size(colourHist,2), numBlobs * 2 - 1);
blobStruct.textureHist = zeros(size(textureHist,2), numBlobs * 2 - 1);
blobStruct.size = zeros(numBlobs * 2 -1, 1);
blobStruct.boxes = zeros(numBlobs * 2 - 1, 4);

% Insert calculated histograms, sizes, and boxes
blobStruct.colourHist(:,1:numBlobs) = colourHist';
blobStruct.textureHist(:,1:numBlobs) = textureHist';
blobStruct.size(1:numBlobs) = blobSizes ./ 3;
blobStruct.boxes(1:numBlobs,:) = blobBoxes;

blobStruct.imSize = size(colourIm,1) * size(colourIm,2);

%%% If you want to use original blobs in similarity functions, uncomment
%%% these lines.
% blobStruct.blobs = cell(numBlobs * 2 - 1, 1);
% initialBlobs = SegmentIndices2Blobs(blobIndIm, blobBoxes);
% blobStruct.blobs(1:numBlobs) = initialBlobs;


% Loop over all merging strategies. Perform them one by one.
boxes = cell(1, length(functionHandles)+1);
priority = cell(1, length(functionHandles) + 1);
hierarchy = cell(1, length(functionHandles));
for i=1:length(functionHandles)
    [boxes{i}, hierarchy{i}, blobStructT, mergeThreshold] = ...
        BlobStruct2HierarchicalGrouping(blobStruct, neighbours, numBlobs, functionHandles{i});
    boxes{i} = boxes{i}(numBlobs+1:end,:);
    priority{i} = (size(boxes{i}, 1):-1:1)';
end

% Also save the initial boxes
i = i+1;
boxes{i} = blobBoxes;
priority{i} = ones(size(boxes{i}, 1), 1) * (size(boxes{1}, 1)+1);

% Concatenate boxes and priorities resulting from the different merging
% strategies
boxes = cat(1, boxes{:});
priority = cat(1, priority{:});
[priority, ids] = sort(priority, 'ascend');
boxes = boxes(ids,:);



