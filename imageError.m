function E = imageError(img,enc,filters,method,params)
% IMAGEERROR Computes local reconstructions g at every pixel in the image
%   and outputs the reconstruction error. 
% 
%   E = IMAGEERROR(img,enc,filters,method) where img is an input image and
%   enc is its encoding, computes the local reconstruction error for all 
%   patches of the shapes defined in the filters cell array. Both img and 
%   enc are HxWxCxR arrays and IMAGEERROR assumes that reconstructing a 
%   patch locally is equivalent to "expanding" a single color (RGB or LAB)
%   triplet across the whole disk area. 
% 
%   E = IMAGEERROR(img,enc,[],method) where method is one of the supported
%   histogram distance metrics, assumes that both img and enc are HxWxCxBxR
%   histograms. 
% 
%   E = IMAGEERROR(img,enc,filters,'smirnov',params) assumes img is a HxWxCxR
%   input image and enc is its HxWxCxBxR histogram encoding. In this case,
%   IMAGEERROR computes local reconstructions using the inverse transform
%   sampling or smirnov method, for creating random samples from a
%   histogram distribution, and then computes the reconstruction error
%   using the method specified in params.
% 
%   The supported methods for computing the reconstruction error are:
% 
%   {(r)se}:  (root) squared error E = sqrt(sum((y-x).^2))
%   (r)mse :  root mean squared error E = sqrt(mean((y-x).^2))
%   n(r)mse:  normalized rms E = sqrt(sum((y-x).^2) / sum(x.^2))
%   dssim:    structural dissimilarity E = (1-ssim(y,x))/2, where ssim is 
%             the structural similarity index. This should be used in the
%             RGB color space.
%   hist-chi2: chi-squared histogram distance (assumes inputs img and enc
%              are histograms of dimensions HxWxCxBxR.
%   hist-chi2-kernel: chi-squared histogram distance + kernel density 
%                     estimate.
%   hist-intersection: 1-histogram intersection distance (assumes inputs 
%              img and enc are histograms of dimensions HxWxCxBxR.
%   
%   See also: immse, ssim, patchError, compressHistogram
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016


switch method
    % ---------------------------------------------------------------------
    % Inverse tranform sampling or Smirnov transform
    % ---------------------------------------------------------------------
    case {'smirnov','inverse'}
        error('Under construction. Not ready yet.')
        [H,W,C,B,R] = size(enc);
        binCenters = (0:(B-1))/B;
        % Compute filter areas at all scales
        areas = zeros(numel(filters));
        for r=1:R, areas(r) = nnz(filters{r}); end
        % Compute cumulative distribution functions
        CDF = cumsum(enc,4);
        for r=1:R
            % Create random numbers in (0,1). Single precision to save RAM. 
            % Maybe preallocate to save time.
            u = rand(H,W,C,1,areas(r),'single');
            % Hijack max function. Its input is logical, so it will return
            % the index of the first element across the selected dimension
            % that is nonzero.
            [~,u] = max(bsxfun(@gt,u,CDF),[],4);
            u = reshape(binCenters(u), H,W,C,[]);
        end
    % ---------------------------------------------------------------------
    % Histogram comparison methods
    % ---------------------------------------------------------------------
    case {'hist-intersection','hist-chi2','hist-chi2-kernel'}
        if strcmp(method,'hist-intersection')
            E = sum(min(img,enc),4); % histogram intersection
        elseif strcmp(method,'hist-chi2')
            E = 0.5*sum((img-enc).^2 ./ (img + enc + eps), 4);
        else
            [H,W,C,B,R] = size(img);
            % Compute color bin weighted distance using gaussian kernel
            binCenters = ((1:B)-0.5)/B;
            [x,y] = meshgrid(binCenters,binCenters);
            % Compute distance at each channel and scale
            elab = zeros(H*W*3,R);
            imgc = reshape(img(:,:,1:3,:,:),H*W*3,B,R);
            encc = reshape(enc(:,:,1:3,:,:),H*W*3,B,R);
            for r=1:R
                binCenterDistance = 1-exp(-(x-y).^2./(2*r.^2)); % BxB
                dabs = abs(imgc(:,:,r)-encc(:,:,r)); % H*W*C x B
                elab(:,r) = sum((dabs*binCenterDistance) .* dabs, 2);
            end
            elab = reshape(elab,H,W,3,1,R);
            % If we use texture, it's the 4th channel.
            % For texture we compute the chi-squared distance as usual.
            if C > 3
                imgt = img(:,:,4,:,:); enct = enc(:,:,4,:,:);
                etexture = 0.5*sum((imgt-enct).^2 ./ (imgt + enct + eps), 4);
            else
                etexture = [];
            end
            E = cat(3,elab,etexture);
        end
        E(:,:,2,:,:) = E(:,:,2,:,:) + E(:,:,3,:,:); % merge color channels
        E(:,:,3,:,:) = [];  % remove redundant color channel
        E = squeeze(E);
    % ---------------------------------------------------------------------        
    % Pixels-wise squared and root squared metrics
    % ---------------------------------------------------------------------
    case {'se','mse','nmse','rse','rmse','nrmse'}
        [H,W,C,R] = size(enc);
        E = zeros(H,W,R);        
        w = ones(1,C)/C; % be default all channels have equal weights
        if nargin == 5 && ~isempty(params)
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
            % As a result, in the end we get Sum(g-I)^2 = Sum(I.^2) - A2*g.^2.
            % The following quantities are sums withing the disk region D.
            for c=1:C
                I = img(:,:,c);
                g = enc(:,:,c,r);
                sumI2 = conv2(I.^2, D, 'same');
                numer(:,:,c) = sumI2 - A * g.^2;                
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
    % ---------------------------------------------------------------------        
    % Structural similarity based metric
    % ---------------------------------------------------------------------        
    case 'dssim'
        [H,W,C,R] = size(enc);
        E = zeros(H,W,R);        
        k1 = 0.01; k2 = 0.03; L = 1; % default constant values
        if nargin == 5 && ~isempty(params)
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
                I   = img(:,:,c);
                g   = enc(:,:,c,r);
                mi  = conv2(I, D/A, 'same');
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