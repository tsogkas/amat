function fh = drawHistGradOnFigureInteractive(imgRGB,filters)
% DRAWHISTGRADONFIGUREINTERACTIVE Draw histogram gradient responses on the 
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
r  = 5;  % default radius
dr = 1;  % default radius difference when comparing histograms
B  = 32; % used for histogram encodings
R  = 25; % default range of scales
channelType   = {'luminance','colora','colorb','texture'};
distanceType  = {'chi2','interection','maxabs'};
distanceIndex = find(strcmp(distanceType, 'chi2'));
channelIndex  = find(strcmp(channelType, 'luminance'));

% Plot figure and set callbacks
fh = figure; subplot(121); imshow(imgRGB);
set(fh, 'WindowButtonDownFcn',   @changeMethod);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
% set(fh, 'WindowButtonMotionFcn', @changePoint);


% Compute Lab color transformation and histograms with default parameters
[H,W,C] = size(imgRGB);
imgRGB  = reshape(imgRGB, [], C);
imgLab  = rgb2labNormalized(imgRGB);
h       = computeHistograms();

% Construct filters if needed
if nargin < 2, 
    filters = cell(R,1); for i=1:R, filters{i} = disk(i); end; 
end




    function drawHistogramGradients(fh)
        % Get point coordinates and check for validity
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        if x < 1 || x > W || y < 1 || y > H
            title('You are outside of the figure'); drawnow; return
        end
        if x-r < 1 || x+r > W || y-r < 1 || y+r > H
            title('Inner disk crosses the image boundary'); drawnow; return
        end
        if x-r-dr < 1 || x+r+dr > W || y-r-dr < 1 || y+r+dr > H
            title('Outer disk crosses the image boundary'); drawnow; return
        end
        
        % Compute distance
        c = channelIndex;
        switch distanceType
            case 'intersection'
                d = 1-sum(min(h(:,:,c,:,r+dr),h(:,:,c,:,r)),4);
            case 'chi2'
                d = 0.5*sum(((h(:,:,c,:,r+dr)-h(:,:,c,:,r)).^2) ./ ...
                    (h(:,:,c,:,r+dr)+h(:,:,c,:,r)+eps), 4);
            case 'maxabs'
                d = max(abs(h(:,:,c,:,r+dr)-h(:,:,c,:,r)),[],4);
            otherwise, error('Distance type is not supported')
        end

        % Disable annoying docking error that clutters the command line
        if strcmp(fh.WindowStyle, 'docked')
            warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
        end
        % Display image and then restore original patch
        subplot(122); imagesc(d); axis off image;
        title(sprintf('Point (%d,%d), r=%d, hist (%d bins), %s: %.4f',...
                x,y,r,B,errorType,err));
        drawnow;
    end


    function h = computeHistograms()
        % Compute histograms
        f = cat(3, binImage(imgLab,B),textonMap(imgRGB, B));
        h = zeros(H,W,size(f,3),B,R);
        for c=1:size(f,3)
            imgc = f(:,:,c);
            for b=1:B
                imgcb = double(imgc == b);
                for s=1:R
                    h(:,:,c,b,s) = conv2(imgcb, ...
                        double(filters{s})/nnz(filters{s}), 'same');
                end
            end
        end
    end

    function changeRadius(fh,callbackData)
        r = min(min(H,W)/2, max(1, r + callbackData.VerticalScrollCount));
        drawHistogramGradients(fh);
    end
    
    function changeHistGradParams(fh)
        validInput = false;
        validNumBins = false;
        validDr = false;
        prompt = {'Enter number of bins:','Enter dr:'};
        dlgTitle = 'Input';
        while ~validInput
            answer = inputdlg(prompt,dlgTitle);
            if isempty(answer)
                validInput = true; % keep previous #bins and dr
            else
                % Check if answer for #bins is valid
                if isempty(answer{1})
                    validNumBins = true;
                else
                    answer1 = str2double(answer{1});
                    if answer1 <= 0 || isnan(answer1) || isinf(answer1)
                        dlgTitle = 'Invalid input! #bins must be a positive scalar.';
                    else
                        B = answer1;
                        validNumBins = true;
                    end
                end
                % Check if answer for dr is valid
                if isempty(answer{2})
                    validDr = true;
                else
                    answer2 = str2double(answer{2});
                    if answer2 <= 0 || isnan(answer2) || isinf(answer2)
                        dlgTitle = 'Invalid input! dr must be a positive scalar.';
                    else
                        dr = answer2;
                        validDr = true;
                    end
                end
                validInput = validNumBins && validDr;
            end
        end
        drawHistogramGradients(fh);
    end

    function changeMethod(fh,~)
        if strcmp(fh.SelectionType, 'normal')  % change channel
            channelIndex = max(1,mod(channelIndex + 1, numel(channelType)+1));
        elseif strcmp(fh.SelectionType, 'alt') % change distance type
            distanceIndex = max(1,mod(distanceIndex + 1, numel(distanceType)+1));
        elseif strcmp(fh.SelectionType, 'extend') % change #bins
            changeHistGradParams(fh);
            h = computeHistograms(); % must recompute the histograms
        end
        drawHistogramGradients(fh)
    end

    function changePoint(fh,~)
        drawDisk(fh);
    end

end