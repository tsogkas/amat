function [rec, comp] = seg2reconstruction(img,seg)
img       = im2double(img);
numSegs   = size(seg,3);
nnzPixels = zeros(numSegs,1);
rec       = zeros([size(img), numSegs]);
for i=1:numSegs
    segi = seg(:,:,i);
    labels = unique(segi);
    numSegments = numel(labels);
    for s=1:numSegments
        segment = segi == labels(s);
        area = nnz(segment);
        boundaries = bwperim(segment);
        nnzPixels(i) = nnzPixels(i) + nnz(boundaries);
        meanColor = sum(sum(bsxfun(@times, img, double(segment))))/area;
        meanSegment = bsxfun(@times, double(repmat(segment,1,1,3)), meanColor);
        rec(:,:,:,i) = rec(:,:,:,i) + meanSegment; 
    end
end
% Select segmentation with maximum compression and return respective rec
comp = min(nnzPixels / (size(img,1)*size(img,2)));
