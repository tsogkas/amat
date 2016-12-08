function [xc,yc,rc] = findContainedDisks(R)
% Find center coordinates and respective radii of all disks that are FULLY
% contained in any given disk of radius R. Note that this function returns
% RELATIVE coordinates with respect to the disk's center.

% Center coordinates
[x,y] = meshgrid(-R:R,-R:R);
x2 = x.^2; y2 = y.^2;
xc = []; yc = []; rc = [];
for r=1:R
   dc = x2 + y2 <= (R-r).^2;
   xc = [xc; x(dc)];
   yc = [yc; y(dc)];
   rc = [rc; r*ones(length(xc),1)];
end
