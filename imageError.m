function E = imageError(img,enc,filters,method)
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
%   {rse}:  root squared error E = sqrt(sum((y-x).^2))
%   rmse :  root mean squared error E = sqrt(mean((y-x).^2))
%   nrmse:  normalized rms E = sqrt(sum((y-x).^2) / sum(x.^2))
%   dssim:  structural dissimilarity E = (1-ssim(y,x))/2, where ssim is the
%           structural similarity index.
%   
%   See also: immse, ssim
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016


[H,W,numChannels,numScales] = size(enc);
E = zeros(H,W,numScales);
switch method
    case {'rse','rmse','nrmse'}
        for r=1:numScales
            if strcmp(method,'nrmse')
                denom = zeros(H,W);
            end
            numer = zeros(H,W);
            D = double(filters{r});
            A = nnz(D);
            % Compute the terms needed to compute the (relative) NRMS error:
            % Sum ||g(x)-I(x)||^2 / Sum||I(x)||^2. We can be more efficient by
            % expanding the identity: abs(g-I)^2 = (g-I)^2 = g^2 + I^2 - 2*g*I
            % Also note that since we expand a single scalar value across the
            % disk area, Sum(g^2) = A*g^2, where A is the area of the disk.
            % The following quantities are sums withing the disk region D.
            for c=1:numChannels
                I = img(:,:,c);
                g = enc(:,:,c,r);
                sumI2 = conv2(I.^2, D, 'same');
                sumI  = conv2(I,    D, 'same');
                sumg2 = A * g.^2;
                sumgI = g .* sumI;
                numer = numer + sumg2 + sumI2 - 2*sumgI;
                if strcmp(method,'nrmse')
                    denom = denom + sumI2;
                end
            end
            if strcmp(method,'rse')            
                E(:,:,r) = max(0,numer);
            elseif strcmp(method,'rmse')
                E(:,:,r) = max(0, min(1, numer ./ (3*A)));
            elseif strcmp(method,'nrmse')
                E(:,:,r) = max(0, min(1, numer ./ denom));
            end
        end        
        E = sqrt(E); 
    case 'dssim'
        % default constant values (wikipedia)
        k1 = 0.01; k2 = 0.03; L  = 1;
        c1 = (k1*L)^2; c2 = (k2*L)^2;
        for r=1:numScales
            D = double(filters{r});
            A = nnz(D);
            ssim = zeros(H,W);
            % Channel-wise implementation of ssim and average across
            % channels. This is a simplified approach of the matlab
            % implementation, that avoids gaussian weighting for
            % efficiency.
            for c=1:numChannels
                I  = img(:,:,c);
                g  = enc(:,:,c,r);
                mi = conv2(I, D/A, 'same');
                si = conv2((I-mi).^2, D/A, 'same');
                mg = g; % since g is constant within the disk area
                sg = 0;
                sig= 0; 
                ssim = ssim + ((2 .* mi .* mg) .* (2 .* sig + c2)) ./ ...
                              ((mg.^2 + mi.^2 + c1).* (si.^2 + sg.^2 + c2));
            end
            E(:,:,r) = (1-ssim/numChannels)/2; % DSSIM
        end
        E = max(-1, min(1,E));
    otherwise, error('Error type not supported')
end