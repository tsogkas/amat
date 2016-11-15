function e = patchError(img,center,radius,encodeMethod,encodeParam)
% TODO: add errorMethod as additional input argument
if nargin < 4, encodeMethod = 'average'; end

p = diskPatch(img,center,radius);
f = patchEncoding(p,encodeMethod,encodeParam);
nomin = sum(sum(bsxfun(@minus,p,f).^2));
denom = sum(sum(p.^2));
e = sqrt(nomin ./ denom);
