function patch = extractImageDiskPatch(img,center,radius)
% img: input image
% center: [row,col] coordinates of the disk center
% radius: disk radius (scalar)
[rows,cols,channels] = size(img);
[x,y] = meshgrid(1:rows,1:cols);
mask = (x - center(2)).^2 + (y - center(1)).^2 <= radius^2;
img = reshape(img,rows*cols,channels);
patch = img(mask(:), :);
