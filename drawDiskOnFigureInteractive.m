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
methods.encoding = {'average','hist-smirnov','hist-expectation','hist-mode'};
encodingType  = 'average';
errorType = 'se';
errorCounter = find(strcmp(methods.error, errorType));
encodingCounter = find(strcmp(methods.encoding, encodingType));

% Plot figure and set callbacks
fh = figure(1); imshow(imgRGB);
set(fh, 'WindowButtonMotionFcn', @changePoint);
set(fh, 'WindowButtonDownFcn',   @changeMethod);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
[H,W,C] = size(imgRGB);
tmap    = textonMap(imgRGB, numBins); 
tmap    = reshape(tmap,H*W,[]);
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
        if strcmp(encodingType,'average')
            encPatch = patchEncoding(imgPatch,'average');
        else
            encPatch = patchEncoding(binImage(imgPatch,numBins),'hist',numBins);
        end
        switch encodingType
            case 'average'
                err = patchError(imgPatch,encPatch,errorType);
            case 'hist-smirnov'
                % for patch encoding using the hist method, we will try to
                % decode the patch with the Smirnov transform.
                hists = encPatch;
                histt = patchEncoding(tmap(D,:),'hist',numBins);
                encPatch = histinv(encPatch,size(imgPatch,1))';
                err = patchError(imgPatch,encPatch,errorType);
            case 'hist-mode'
                hists = encPatch;
                histt = patchEncoding(tmap(D,:),'hist',numBins);
                encPatch = reshape(compressHistogram(hists,'mode',2),1,[]);
                err = patchError(imgPatch,encPatch,errorType);
            case 'hist-expectation'
                hists = encPatch;
                histt = patchEncoding(tmap(D,:),'hist',numBins);
                encPatch = reshape(compressHistogram(hists,'expectation',2),1,[]);
                err = patchError(imgPatch,encPatch,errorType);
            otherwise, error('Encoding type not supported')
        end
        
        
        % Replace pixels in the input image
        originalPatch = imgRGB(D,:);
        if ~strcmp(errorType, 'dssim')
            encPatch = labNormalized2rgb(encPatch);
        end
        if isvector(encPatch) && size(encPatch,2) == 3
            imgRGB(D,:) = repmat(encPatch, [nnz(D), 1]);
        elseif ismatrix(encPatch)
            imgRGB(D,:) = encPatch;
        else error('Something went wrong with encPatch')
        end
        
        % Disable annoying docking error that clutters the command line
        if strcmp(fh.WindowStyle, 'docked')
            warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
        end
        % Display image and then restore original patch
        figure(1); imshow(reshape(imgRGB,H,W,C)); imgRGB(D,:) = originalPatch; drawnow;
        if strcmp(encodingType, 'average')
            title(sprintf('Point (%d,%d), r=%d, average, %s: %.4f',...
                x,y,r,errorType,err));
        else
            title(sprintf('Point (%d,%d), r=%d, %s (%d bins), %s: %.4f',...
                x,y,r,encodingType,numBins,errorType,err));
            figure(2);
            if C == 1 
                subplot(121); bar(hists(1,:)); title('luminance');
                subplot(122); bar(histt); title('texture')
            else
                subplot(221); bar(hists(1,:)); title('luminance');
                subplot(222); bar(hists(2,:)); title('colora');
                subplot(223); bar(hists(3,:)); title('colorb');
                subplot(224); bar(histt); title('texture');
            end
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
end