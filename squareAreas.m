function out = squareAreas(H,W,r,deg)
%SQUAREAREA
%   returns an HxW matrix that represents how many pixels are covered by a
%   square with radius r
%   rotation deg
if nargin < 4
    deg = 0;
end

rotatedFilter = imrotate(ones(2*r+1), deg);
rotatedFilter = padarray(rotatedFilter, double(mod(size(rotatedFilter), 2) == 0), 'pre');

out = conv2(ones(H,W), rotatedFilter, 'same');
end

