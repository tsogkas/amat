function E = imageError(img,enc,filters,method,params)
% IMAGEERROR Computes local reconstructions g at every pixel in the image
%   and outputs the reconstruction error. 
% 
%   E = IMAGEERROR(img,enc,filters,method) where img is an input image and
%   enc is its encoding, computes the local reconstruction error for all 
%   patches of the shapes defined in the filters cell array. The function
%   assumes that the reconstruction is equivalent to "expanding" a single
%   RGB triplet across the whole disk area. 
% 
%   The supported methods for computing the reconstruction error are:
% 
%   {(r)se}:  (root) squared error E = sqrt(sum((y-x).^2))
%   (r)mse :  root mean squared error E = sqrt(mean((y-x).^2))
%   n(r)mse:  normalized rms E = sqrt(sum((y-x).^2) / sum(x.^2))
%   dssim:    structural dissimilarity E = (1-ssim(y,x))/2, where ssim is 
%             the structural similarity index. This should be used in the
%             RGB color space.
%   
%   See also: immse, ssim
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016


[H,W,C,R] = size(enc);
E = zeros(H,W,R);
switch method
    case {'se','mse','nmse','rse','rmse','nrmse'}
        w = [1/3 1/3 1/3]; % default channel weights
        if nargin == 5
            assert(sum(params(:)) == 1, 'Channel weights should sum to 1')
            w = params;
        end
        numer = zeros(H,W,C);
        denom = zeros(H,W,C);
        for r=1:R
            numer(:) = 0; denom(:) = 0;
            D = double(filters{r});
            A = nnz(D);
            % Compute the per-channel square error:
            % Sum ||g(x)-I(x)||^2 / Sum||I(x)||^2. We can be more efficient by
            % expanding the identity: abs(g-I)^2 = (g-I)^2 = g^2 + I^2 - 2*g*I
            % Also note that since we expand a single scalar value across the
            % disk area, Sum(g^2) = A*g^2, where A is the area of the disk.
            % The following quantities are sums withing the disk region D.
            for c=1:C
                I = img(:,:,c);
                g = enc(:,:,c,r);
                sumI2 = conv2(I.^2, D, 'same');
                sumI  = conv2(I,    D, 'same');
                sumg2 = A * g.^2;
                sumgI = g .* sumI;
                numer(:,:,c) = sumg2 + sumI2 - 2*sumgI;                
                if strcmp(method,'nrmse') || strcmp(method,'nmse')
                    denom(:,:,c) = sumI2;
                end
            end
            % Normalize across each channel
            if strcmp(method,'rmse') || strcmp(method,'mse')
                numer = numer ./ A;                 
            elseif strcmp(method,'nrmse') || strcmp(method,'nmse')
                numer = numer ./ denom;
            end
            % Combine channels with different weights
            E(:,:,r) = reshape(reshape(numer,[],C)*w', H,W);
        end
        % Make positive and take square root if needed
        E = max(0,E);
        if any(strcmp(method,{'rse','rmse','nrmse'}))
            E = sqrt(E); 
        end
    case 'dssim'
        k1 = 0.01; k2 = 0.03; L = 1; % default constant values
        if nargin == 5
            k1 = params(1); k2 = params(2); L  = params(3);
        end
        c1 = (k1*L)^2; c2 = (k2*L)^2;
        for r=1:R
            D = double(filters{r});
            A = nnz(D);
            ssim = zeros(H,W);
            % Channel-wise implementation of ssim and average across
            % channels. This is a simplified approach of the matlab
            % implementation, that avoids gaussian weighting for
            % efficiency.
            for c=1:C
                I  = img(:,:,c);
                g  = enc(:,:,c,r);
                mi = conv2(I, D/A, 'same');
                si2 = conv2((I-mi).^2, D/(A-1), 'same'); % unbiased estimate
                mg  = g; % since g is constant within the disk area
                sg2 = 0;
                sig = 0; 
                ssim = ssim + ((2 .* mi .* mg) .* (2 .* sig + c2)) ./ ...
                              ((mg.^2 + mi.^2 + c1).* (si2 + sg2 + c2));
            end
            E(:,:,r) = (1-ssim/C)/2; % DSSIM
        end
        E = max(-1, min(1,E));
    otherwise, error('Error type not supported')
end