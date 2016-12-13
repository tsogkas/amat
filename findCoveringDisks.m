function [xc,yc,rc] = findCoveringDisks(p,R,sz)
% Find center coordinates and respective radii of all disks that cover
% point p = (xp,yp). The maximum radius of a covering disk is R and the
% (optional) size of the image, sz.

% Center coordinates
xp = c(1); yp = c(2);
% Grid for a disk of radius R
[x,y] = meshgrid(-R:R,-R:R);
% Array centers of contained disks at each scale
dc = bsxfun(@le,x.^2 + y.^2,reshape((1:R).^2, 1,1,[]));
[yc,xc,rc] = ind2sub([2*R+1,2*R+1,R],find(dc));
yc = yc - R - 1 + yp; xc = xc - R - 1 + xp;
% If image size (sz) is given, remove invalid points that cross the boundary
if nargin == 3
    outOfLimits = xc < 1 | yc < 1 | xc > W | yc > H;
    xc(outOfLimits) = []; yc(outOfLimits) = []; rc(outOfLimits) = [];
end