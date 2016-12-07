function [xc,yc,rc] = findContainedDisks(c,r)

% Center coordinates
[x,y] = meshgrid(-r:r,-r:r);
x2 = x.^2; y2 = y.^2;
xc = []; yc = []; rc = [];
for rr=1:r
   dc = x2 + y2 <= (r-rr).^2;
   xc = [xc; x(dc)];
   yc = [yc; y(dc)];
   rc = [rc; rr*ones(length(xc),1)];
end
xc = xc + c(1);
yc = yc + c(2);
