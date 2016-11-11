function f = computeDiskEncodings(img,filters)
% In the simplest case, an encoding at point (x,y,r) is the average of the
% values within a disk region of radius r, centered at point (x,y).

[H,W,numChannels] = size(img);
numScales = numel(filters);
circleSum = zeros(H,W,numChannels,numScales);
for c=1:numChannels
    for s=1:numScales
        circleSum(:,:,c,s) = conv2(img(:,:,c), double(filters{s}), 'same');
    end
end
% The disk area is computed as the cumulative sum of perimeters of circles
% of increasing  radii.
circlePerimeter = cellfun(@nnz,filters);
diskArea = cumsum(circlePerimeter);

% The sum within a disk area is computed as the cumulative sum within
% circles of increasing radii.
f = cumsum(circleSum,4); % (x,y,k,r)
f = bsxfun(@rdivide,f,reshape(diskArea,1,1,1,[]));
assert(isinrange(f,[0,1]))