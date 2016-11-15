function p = rgb2disk(rgb,r)
pside = 2*r+1;
p = zeros(pside,pside,3);
[x,y] = meshgrid(-r:r);
p = reshape(p,[],3);
mask = x.^2 + y.^2 <= r^2;
p(mask,:) = repmat(rgb, [nnz(mask),1]);
p = reshape(p,pside,pside,3);