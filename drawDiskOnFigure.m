function h = drawDiskOnFigure(h,~,img,enc,r,method)
% DRAWDISKONFIGURE Interactive drawing of reconstructed disks on the input
%   figure and display of useful info. This function is used as a callback
%   for a certain window action (mouse motion, mouse click) on the figure,
%   and takes as input the following arguments:
% 
%   img:    the input image on which we draw the reconstructed disks.
%   enc:    array with encodings at all pixels and radii.
%   r:      radius of the disk to be drawn.
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

% Get point coordinates and check for validity
x = round(h.CurrentAxes.CurrentPoint(1,1));
y = round(h.CurrentAxes.CurrentPoint(1,2));
[H,W,C] = size(img);
if x < 1 || x > W || y < 1 || y > H
    disp('You clicked outside of the figure')
    return
end

[xx,yy] = meshgrid(1:W,1:H);
mask = (xx-x).^2 + (yy-y).^2 <= r^2;
img = reshape(img, [], C);
img(mask,:) = repmat(reshape(enc(y,x,:,r), [1 1 3]), [nnz(mask),1]);
img = reshape(img,H,W,C);
imshow(img); title(['Point (' num2str(x) ',' num2str(y) ')'])




