function e = patchError(img,center,radius,encodeMethod,encodeParam)
% TODO: update error methods
if nargin < 4, encodeMethod = 'average'; end

p = diskPatch(img,center,radius);
switch encodeMethod 
    case 'average'
        f = patchEncoding(p,encodeMethod);
    case 'hist'
        f = patchEncoding(p,encodeMethod, encodeParam);
    otherwise, error('Encoding method not supported')
end
numer = sum(sum(bsxfun(@minus,p,f).^2));
denom = sum(sum(p.^2));
e = sqrt(numer ./ denom);
