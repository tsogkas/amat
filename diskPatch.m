function p = diskPatch(img,center,radius,returnBox)
% DISKPATCH returns the pixels corresponding to a disk area of a
%   given center and radius in the image.
% 
%   p = DISKPATCH(img,center,radius) where center is a [row,col]
%   vector and radius is a scalar, returns the corresponding disk area of
%   the image. If the input is grayscale, then the output is 
%   numel(diskArea) x 1, otherwise it is numel(diskArea) x 3.
% 
%   p = DISKPATCH(img,center,radius,returnBox) where returnBox is a
%   logical variable (true/false) controls if DISKPATCH will return
%   the rectangular box that encloses the disk, instead of returning just
%   the pixels inside the disk. Default is false.
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 4, returnBox = false; end

[H,W,C] = size(img);
x = center(1); y = center(2); r = radius;
% Check boundaries and radius
if x < 1 || x > W || y < 1 || y > H
    error('Disk center falls out of the image!')
end
if x-r < 1 || x+r > W || y-r < 1 || y+r > H
    error('Disk crosses the image boundary!')
end

% Whatever is outside the image domain is zero.
p = img(y-r:y+r,x-r:x+r,:);
if ~returnBox
    [xx,yy] = meshgrid(-r:r);
    p = reshape(p,[],C);
    p = p(xx.^2 + yy.^2 <= r^2, :);
end


% Old version
% % Limits of the patch within the image domain
% xmin = max(x-radius,1);
% ymin = max(y-radius,1);
% xmax = min(x+radius,W);
% ymax = min(y+radius,H);
% % Whatever is outside the image domain is zero.
% if x-radius >= 1, pxs = 1; else pxs = 2-x+radius; end
% if y-radius >= 1, pys = 1; else pys = 2-y+radius; end
% if x+radius <= W, pxe = pside; else pxe = pside + W - (x+radius); end
% if y+radius <= H, pye = pside; else pye = ps1ide + H - (y+radius); end
% p(pys:pye,pxs:pxe,:)  = img(ymin:ymax,xmin:xmax,:);
% if ~returnBox
%     [xx,yy] = meshgrid(-radius:radius);
%     p = reshape(p,[],C);
%     p = p(xx.^2 + yy.^2 <= radius^2, :);
% end