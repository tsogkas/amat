function out = computeSquareSumsConv(inputIntegral, r)
%COMPUTESQUARESUMSCONV
%   inputIntegral: the integral image of an input I, geneated by calling
%                  integralImage(I)
%   r: radius
%   returns the square sums of the square at each pixel with radius r

fil = zeros(2*r+2, 1);
fil(1) = 1;
fil(end) = -1;

[H,W,C] = size(inputIntegral);
out = zeros(H,W,C);

for c=1:C
    % separable filter
    out(:,:,c) = conv2(conv2(inputIntegral(:,:,c), fil, 'same'), fil', 'same');
end

% make the output the same as input image
out = out(1:end-1, 1:end-1, :);
end

