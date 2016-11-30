function E = maximalityError(enc,dr,method)
% takes as input encoding enc HxWxCxBxR and returns the error of using disk
% of radius r instead of r+dr at each point. This error term should be
% independend of the image reconstruction error at the respective points.

if nargin < 2, dr = 1; end
if nargin < 3, method = 'chi2'; end

[H,W,C,B,R] = size(enc);
E = ones(H,W,C,R); % start by setting all the errors to the max value
for r=1:R
    E(:,:,:,r) = histogramDistance(enc(:,:,:,:,r+dr), enc(:,:,:,:,r),method);
end





