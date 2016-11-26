function e = patchError(imgPatch,encPatch,errorType,errorParams)
% PATCHERROR Computes the error between an encoded patch and the respective
%   patch from the original image.
% 
%   e = PATCHERROR(imgPatch,encPatch,errorType) where imgPatch represents
%   an image patch and encPatch represents its encoding, computes the error
%   of using the metric determined by errorType. PATCHEERROR assumes that
%   imgPatch is a NxC matrix, where N is the number of pixels considered in
%   the patch (which is not necessarily square) and C is the number of
%   channels at each pixel. encPatch is a 1xC vector representing the mean
%   value encoding of the patch, and for the reconstruction of the patch,
%   we consider the repetition of this 1xC tuple across the whole patch.
% 
%   If errorType is one of the histogram distance methods, imgPatch and
%   encPatch are assumed to be the BxC histogram representations. imgPatch
%   is the full histogram representation of the original image patch, while
%   encPatch is its possibly "compressed" version.
% 
%   The available error metrics are the following:
% 
%   {(r)se}:  (root) squared error E = sqrt(sum((y-x).^2))
%   (r)mse :  root mean squared error E = sqrt(mean((y-x).^2))
%   n(r)mse:  normalized rms E = sqrt(sum((y-x).^2) / sum(x.^2))
%   dssim}    computes the local structural dissimilarity (1-ssim(y,x))/2. 
%             This should be used in the RGB color space.
%   hist-chi2: chi-squared histogram distance (assumes inputs img and enc
%              are histograms of dimensions HxWxCxBxR.
%   hist-chi2-kernel: chi-squared histogram distance + kernel density 
%                     estimate.
%   hist-intersection: 1-histogram intersection distance (assumes inputs 
%              img and enc are histograms of dimensions BxC.
% 
%   See also: immse, ssim, imageError, compressHistogram
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, errorType = 'se'; end

switch errorType
    % ---------------------------------------------------------------------
    % Histogram comparison methods
    % ---------------------------------------------------------------------
    case {'hist-intersection','hist-chi2','hist-chi2-kernel'}
        if strcmp(method,'hist-intersection')
            e = sum(min(imgPatch,encPatch)); % histogram intersection
        elseif strcmp(method,'hist-chi2')
            e = 0.5*sum((imgPatch-encPatch).^2 ./ (imgPatch+encPatch+eps), 4);
        else
            if nargin == 4 && ~isempty(errorParams)
                r = errorParams{1};
            else
                error('You must provide the kernel scale (sigma)')
            end
            [B,C] = size(imgPatch);
            % Compute color bin weighted distance using gaussian kernel
            binCenters = ((1:B)-0.5)/B;
            [x,y] = meshgrid(binCenters,binCenters);
            % Estimate kernel scale (sigma) from the number of pixels 
            binCenterDistance = 1-exp(-(x-y).^2./(2*r.^2)); % BxB
            dabs = abs(imgPatch(:,1:3)-encPatch(:,1:3)); 
            elab = sum((binCenterDistance*dabs) .* dabs);
            % If we use texture, it's the 4th channel.
            % For texture we compute the chi-squared distance as usual.
            if C > 3
                imgt = imgPatch(:,4); enct = encPatch(:,4);
                etexture = 0.5*sum((imgt-enct).^2 ./ (imgt + enct + eps));
            else
                etexture = [];
            end
            e = [elab,etexture];
        end
        % merge color channels and remove redundant color channel
        e(2) = e(2) + e(3); e(3) = [];         
    % ---------------------------------------------------------------------        
    % Pixels-wise squared and root squared metrics
    % ---------------------------------------------------------------------
    case {'se','mse','nmse','rse','rmse','nrmse'} 
        % errorParams should be a 1x3 vector containing the weights for the
        % three color channels. This is useful when using a color space
        % that is different than RGB (e.g. CIE Lab), where the channels do
        % not have equal discriminative importance.
        [N,C] = size(imgPatch);
        w = ones(1,C)/C; % be default all channels have equal weights
        if nargin == 4 && ~isempty(errorParams)
            assert(sum(errorParams(:))==1, 'Color channel weights should sum to 1');
            w = errorParams; 
        end
        % Compute per-channel squared error
        e = sum(bsxfun(@minus,imgPatch,encPatch).^2);
        % Normalize across each channel
        if strcmp(errorType,'rmse') || strcmp(errorType,'mse')
            e = e ./ N;
        elseif strcmp(errorType,'nrmse') || strcmp(errorType,'nmse')
            e = e ./ sum(imgPatch.^2);
        end
        % Combine channels, make positive and get sqrt if needed
        e = max(0, w * e');
        if any(strcmp({'rse','rmse','nrmse'},errorType))
            e = sqrt(e);
        end
    % ---------------------------------------------------------------------        
    % Structural similarity based metric
    % ---------------------------------------------------------------------        
    case 'dssim' % NOTE: dssim should be used directly on the RGB space
        k1  = 0.01; k2 = 0.03; L = 1; % default values
        if nargin == 4 && ~isempty(errorParams)
            k1 = errorParams(1); k2 = errorParams(2); L  = errorParams(3);
        end
        c1  = (k1*L)^2; c2 = (k2*L)^2;
        % Channel-wise implementation of ssim
        mx  = mean(imgPatch);
        my  = encPatch;
        sx2 = mean(bsxfun(@minus,imgPatch,mx).^2);
        sy2 = 0;
        sxy = 0;
        e = ((2 .* mx .* my) .* (2 .* sxy + c2)) ./ ... % ssim
            ((mx.^2 + my.^2 + c1).* (sx2 + sy2 + c2));
        e = (1-mean(e,2))/2; % mean across channels and dssim
        e = max(-1, min(1,e));
    otherwise, error('Error type not supported')
end
