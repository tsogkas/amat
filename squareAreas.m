function out = squareAreas(H,W,r)
%SQUAREAREA
%   returns an HxW matrix that represents how many pixels are covered by a
%   square with radius r
out = conv2(ones(H,W), ones(2*r+1), 'same');
end

