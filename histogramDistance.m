function d = histogramDistance(h1,h2,method,r)

assert(isequal(size(h1),size(h2)), 'Histogram dimensions missmatch')
if isvector(h1), [~,dim] = max(size(h1)); % works for column and row vectors
elseif ndims(h1) == 5, dim = 4;  % HxWxCxBxR
else dim = ndims(h1); % in all other cases, bins are the last dim
end

switch method
    case 'intersection'
        d = 1-sum(min(h1,h2),dim); % 1-histogram intersection
    case 'chi2'
        d = 0.5*sum((h1-h2).^2 ./ (h1 + h2 + eps), dim);
    case 'chi2-kernel'
        B = size(h1,dim);
        % Compute color bin weighted distance using gaussian kernel
        binCenters = ((1:B)-0.5)/B;
        [x,y] = meshgrid(binCenters,binCenters);
        % Compute distance at a single scale r (sigma of the gaussian)
        h1 = reshape(h1,[],B);
        h2 = reshape(h2,[],B);
        binCenterDistance = 1-exp(-(x-y).^2./(2*r.^2)); % BxB
        dabs = abs(h1-h2); % H*W*C x B
        d = sum((dabs*binCenterDistance) .* dabs, 2);
    otherwise, error('Histogram comparison method not supported')
end