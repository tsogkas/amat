function out = computeSquareSumsConv(inputIntegral, r)
%COMPUTESQUARESUMSCONV
%   inputIntegral: the integral image of an input I, geneated by calling
%                  integralImage(I)
%   r: radius
%   returns the square sums of the square at each pixel with radius r

fil = zeros(2*r+3, 1);
fil(1) = 1;
fil(end-1) = -1;

% pad the matrix to handle the border
paddedArray = padarray(inputIntegral, [r, r], 'replicate', 'post');

[H,W,C] = size(paddedArray);
out = zeros(H,W,C);

for c=1:C
    % separable filter
    out(:,:,c) = conv2(conv2(paddedArray(:,:,c), fil, 'same'), fil', 'same');
end

% make the output the same as input image
out = out(1:end-1-r, 1:end-1-r, :);
end

