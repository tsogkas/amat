function fh = drawDiskOnFigureInteractive(imgRGB)
% DRAWDISKONFIGUREINTERACTIVE Draw reconstructed disk patches on the 
%   original image and display useful information. This tool allows the 
%   user to interactively change various parameters, such as the size of
%   the disk, the type of the error function used, and the method used for
%   summarizing an image patch. The supported functions include:
% 
%   [hover mouse over figure]: change the coordinates of the disk center.
%   [left click]: change the error type used ({'se'},'mse','nmse','rse',
%                 'rmse','nrmse','dssim').
%   [right click]:  change the encoding method ({'average'},'hist').
%   [middle click]: change the number of bins used for the histogram
%                   computations.
%   [scroll wheel]: up/down (decrease/increase) the radius of the disk.
%
%   See also: patchEncoding
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

% Default parameters
r = 5;
numBins = 32; % used for histogram encodings
methods.error    = {'se','mse','nmse','rse','rmse','nrmse','dssim'};
methods.encoding = {'average','hist'};
encodingType  = 'average';
errorType = 'se';
errorCounter = find(strcmp(methods.error, errorType));
encodingCounter = find(strcmp(methods.encoding, encodingType));

% Plot figure and set callbacks
fh = figure; imshow(imgRGB);
set(fh, 'WindowButtonMotionFcn', @changePoint);
set(fh, 'WindowButtonDownFcn',   @changeMethod);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
[H,W,C] = size(imgRGB);
imgRGB  = reshape(imgRGB, [], C);
imgLab  = rgb2labNormalized(imgRGB);
[xx,yy] = meshgrid(1:W,1:H);


    function drawDisk(fh)
        % Get point coordinates and check for validity
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        if x < 1 || x > W || y < 1 || y > H
            title('You are outside of the figure'); drawnow; return
        end
        if x-r < 1 || x+r > W || y-r < 1 || y+r > H
            title('Disk crosses the image boundary'); drawnow; return
        end
        
        % Disk logical mask
        D = (xx-x).^2 + (yy-y).^2 <= r^2;

        % The dssim metric should be used on RGB data
        if strcmp(errorType, 'dssim')
            imgPatch = imgRGB(D,:);
        else
            imgPatch = imgLab(D,:);
        end
        
        % Encode the patch and compute error
        encPatch = patchEncoding(imgPatch,encodingType,numBins);
        err = patchError(imgPatch,encPatch,errorType,[.7 .15 .15]);
        
        % Replace pixels in the input image
        originalPatch = imgRGB(D,:);
        if ~strcmp(errorType, 'dssim')
            encPatch = labNormalized2rgb(encPatch);
        end
        imgRGB(D,:) = repmat(encPatch, [nnz(D), 1]);
        
        % Disable annoying docking error that clutters the command line
        if strcmp(fh.WindowStyle, 'docked')
            warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
        end
        % Display image and then restore original patch
        imshow(reshape(imgRGB,H,W,C)); imgRGB(D,:) = originalPatch;         
        if strcmp(encodingType, 'hist')
            title(sprintf('Point (%d,%d), r=%d, hist (%d bins), %s: %.4f',...
                x,y,r,numBins,errorType,err));
        else
            title(sprintf('Point (%d,%d), r=%d, average, %s: %.4f',...
                x,y,r,errorType,err));
        end
        drawnow;
    end


    function changePoint(fh,~)
        drawDisk(fh);
    end

    function changeRadius(fh,callbackData)
        r = min(min(H,W)/2, max(1, r + callbackData.VerticalScrollCount));
        drawDisk(fh);
    end
    
    function changeNumBins(fh)
        validInput = false;
        dlgTitle = 'Change number of histogram bins';
        while ~validInput
            answer = inputdlg('Enter number of bins:',dlgTitle);
            if isempty(answer)
                validInput = true; % keep previous nBins
            else
                answer = answer{1};
                answer = str2double(answer);
                if isempty(answer) || answer <= 0
                    dlgTitle = 'Invalid input! #bins must be a positive scalar.';
                else
                    numBins = answer;
                    validInput = true;
                end
            end
        end
        drawDisk(fh);
    end

    function changeMethod(fh,~)
        if strcmp(fh.SelectionType, 'normal')
            errorCounter = max(1,mod(errorCounter + 1, numel(methods.error)+1));
            errorType = methods.error{errorCounter};
        elseif strcmp(fh.SelectionType, 'alt')
            encodingCounter = max(1,mod(encodingCounter + 1, numel(methods.encoding)+1));
            encodingType = methods.encoding{encodingCounter};
        elseif strcmp(fh.SelectionType, 'extend') 
            changeNumBins(fh)
        end
        drawDisk(fh)
    end

%     function e = patchError(imgPatch,encPatch)
%         switch errorType
%             % These should be computed in the CIE Lab color space.
%             case {'se','mse','nmse','rse','rmse','nrmse'}
%                 e = sum(sum(bsxfun(@minus,imgPatch,encPatch).^2));
%                 % Normalize
%                 if strcmp(errorType,'rmse') || strcmp(errorType,'mse')
%                     e = e / (C*size(imgPatch,1));
%                 elseif strcmp(errorType,'nrmse') || strcmp(errorType,'nmse')
%                     e = e / sum(imgPatch(:).^2);
%                 end
%                 e = max(0,e);
%                 if ismember(errorType,{'rse','rmse','nrmse'})
%                     e = sqrt(e); 
%                 end                
%             case 'dssim' % NOTE: dssim operates on the RGB color space
%                 % default constant values (wikipedia)
%                 k1 = 0.01; k2 = 0.03; L  = 1;
%                 c1 = (k1*L)^2; c2 = (k2*L)^2;
%                 % Channel-wise implementation of ssim
%                 mx = mean(imgPatch);
%                 sx = mean(bsxfun(@minus,imgPatch,mx).^2);
%                 my = encPatch;
%                 sy = 0;
%                 sxy= 0;
%                 e = ((2 .* mx .* my) .* (2 .* sxy + c2)) ./ ... % ssim
%                     ((mx.^2 + my.^2 + c1).* (sx.^2 + sy.^2 + c2));
%                 e = (1-mean(e,2))/2; % mean across channels and dssim
%                 e = max(-1, min(1,e));
%             otherwise, error('Error type not supported')
%         end
%     end

end