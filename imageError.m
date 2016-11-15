function E = reconstructionError(img,f,filters)
% DECODEIMAGE Computes local reconstructions g at every pixel in the image
%   and outputs the reconstruction error. For now we are using the NRMS
%   error which is not ideal.
% 
%   TODO: add implementation for ssim metric.
% 
%   INPUT:
%   f: HxWxCxR encoding array
%   img: original image
%   method: temporarily used for profiling. either 'conv' or 'im2col'. Need
%   to see which is faster/more efficient
% 
%   OUTPUT:
%   e: HxWxR error array. e(i,j,r) is the reconstruction error of a
%   disk of radius r, placed at pixel (i,j)
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016


[H,W,numChannels,numScales] = size(f);
E = zeros(H,W,numScales);
for r=1:numScales
    nomin = zeros(H,W);
    denom = zeros(H,W);
    D = double(filters{r});
    A = nnz(D);
    for c=1:numChannels
        % Compute the terms needed to compute NRMS error:
        % Sum ||g(x)-I(x)||^2 / Sum||I(x)||^2. We can be more efficient by
        % expanding the identity: abs(g-I)^2 = (g-I)^2 = g^2 + I^2 - 2*g*I
        % Also note that since we expand a single scalar value across the
        % disk area, Sum(g^2) = A*g^2, where A is the area of the disk.
        % The following quantities are sums withing the disk region D.
        I = img(:,:,c);
        g = f(:,:,c,r);
        sumI2 = conv2(I.^2, D, 'same');
        sumI  = conv2(I,    D, 'same');
        sumg2 = A * g.^2;
        sumgI = g .* sumI;
        nomin = nomin + sumg2 + sumI2 - 2*sumgI;
        denom = denom + sumI2;
    end
    E(:,:,r) = nomin ./ denom;
end
E = sqrt(max(0, min(1,E))); 
