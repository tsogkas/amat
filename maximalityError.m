function E = maximalityError(enc,dr)
% takes as input encoding enc HxWxCxBxR and returns the error of using disk
% of radius r instead of r+dr at each point. This error term should be
% independend of the image reconstruction error at the respective points.

if nargin < 2, dr = 1; end
E = imageError(enc(:,:,:,:,1:(R-dr)), enc(:,:,:,:,(dr+1):R));






