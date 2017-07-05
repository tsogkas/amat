function rotatedSqrSum = computeRotatedSquareSums(integralColumn, radius, deg)
%COMPUTEROTATEDSQUARESUMS Summary of this function goes here
%   Detailed explanation goes here
% if the image is I, then integralColumn should be cumsum(I, 1)

[H,W,C] = size(integralColumn);

filter = ones(2*radius+1);

% rotated square
rotatedFilter = imrotate(filter, deg);

% make sure both the width and height of rotatedFilter are odd to make sure
% it has a center
rotatedFilter = padarray(rotatedFilter, double(mod(size(rotatedFilter), 2) == 0), 'pre');
padSize = floor(size(rotatedFilter) / 2);

% pad integralColumn to deal with pixels along the border
cumA = padarray(padarray(integralColumn, padSize, 'pre'), [padSize(1),0], 'post', 'replicate');

simpleRotatedFilter = getSimpleRotatedFilter(rotatedFilter);

rotatedSqrSum = zeros(H,W,C);
for c=1:C
    rotatedSqrSumC = conv2(cumA(:,:,c), simpleRotatedFilter, 'same');
    % make sure the output size is consistent
    rotatedSqrSum(:,:,c) = rotatedSqrSumC(padSize(1):end-padSize(1)-1, padSize(2)+1:end);
end
end
