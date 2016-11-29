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
%   The available error metrics are the following:
% 
%   {(r)se}:  (root) squared error E = sqrt(sum((y-x).^2))
%   (r)mse :  root mean squared error E = sqrt(mean((y-x).^2))
%   n(r)mse:  normalized rms E = sqrt(sum((y-x).^2) / sum(x.^2))
%   dssim}    computes the local structural dissimilarity (1-ssim(y,x))/2. 
%             This should be used in the RGB color space.
% 
%   See also: immse, ssim, imageError
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if nargin < 3, errorType = 'se'; end

switch errorType
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
