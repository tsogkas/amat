function patch = extractDiskPatch(img,center,radius,returnBox)
% EXTRACTDISKPATCH returns the pixels corresponding to a disk area of a
%   given center and radius in the image.
% 
%   p = EXTRACTDISKPATCH(img,center,radius) where center is a [row,col]
%   vector and radius is a scalar, returns the corresponding disk area of
%   the image. If the input is grayscale, then the output is 
%   numel(diskArea) x 1, otherwise it is numel(diskArea) x 3.
% 
%   p = EXTRACTDISKPATCH(img,center,radius,returnBox) where returnBox is a
%   logical variable (true/false) controls if EXTRACTDISKPATCH will return
%   the rectangular box that encloses the disk, instead of returning just
%   the pixels inside the disk.
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 4, returnBox = false; end

[rows,cols,channels] = size(img);
[x,y] = meshgrid(1:rows,1:cols);
mask = (x - center(2)).^2 + (y - center(1)).^2 <= radius^2;
if returnBox
    patch = cropImageBox(img,mask2bbox(mask));
else
    img   = reshape(img,rows*cols,channels);
    patch = img(mask(:), :);
end