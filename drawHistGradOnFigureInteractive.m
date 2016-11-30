function fh = drawHistGradOnFigureInteractive(imgRGB,filters)
% DRAWHISTGRADONFIGUREINTERACTIVE Draw "gradients of histograms" responses 
%   between two homocentric disks at every pixel in the image. This
%   interactive tool allows the user to interactively change various 
%   parameters, such as the size of the inner and outer disks, the type of 
%   the distance metric used, and the number of bins used to compute
%   histograms. The supported functions include:
% 
%   [left click]: change the channel visualized
%                 {'luminance'},'colora','colorb','texture'. For color
%                 perform all computations in the Lab color space and for
%                 textons we used the method and code from [1].
%   [right click]:  change the distance metric used {'chi2'},'intersection'.
%   [middle click]: change the number of bins used for the histogram
%                   computations, and the difference between inner and
%                   outer disk radii, dr.
%   [scroll wheel]: up/down (decrease/increase) the radius of the inner disk.
%
%   [1] Martin, D., Fowlkes, C., and Malik, J.,
%   Learning to detect natural image boundaries using local brightness,
%   color, and texture cues. IEEE Trans. PAMI, 2004.
% 
%   See also: drawDiskOnFgureInteractive, imageEncoding
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

% Default parameters
r  = 1;  % default radius
dr = 1;  % default radius difference when comparing histograms
B  = 32; % used for histogram encodings
R  = 40; % default range of scales
channelType   = {'luminance','color','texture'};
distanceType  = {'chi2','chi2-gaussian','reconstruction','combined'};
distanceIndex = find(strcmp(distanceType, 'chi2'));
channelIndex  = find(strcmp(channelType, 'luminance'));

% Construct filters if needed
if nargin < 2, 
    filters = cell(R,1); for i=1:R, filters{i} = disk(i); end; 
end

% Compute Lab color transformation and histograms with default parameters
[H,W,~] = size(imgRGB);
imgLab  = rgb2labNormalized(imgRGB);
imgBinned = cat(3, binImage(imgLab,B),textonMap(imgRGB, B));
h = imageEncoding(imgBinned,filters,'hist',B);
m = imageEncoding(imgLab,filters,'average');

% Plot figure and set callbacks
fh = figure; subplot(121); imshow(imgRGB); 
title(sprintf('Input image (%dx%d)',H,W));
set(fh, 'WindowButtonDownFcn',   @changeMethod);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
drawHistogramGradients(fh); % first draw

    function drawHistogramGradients(fh)
        % Compute distance
        et = 'mse';
        dr = ceil(r/3);
        if channelIndex == 2, c = [2;3]; else c = channelIndex; end
        switch distanceType{distanceIndex}
            case 'intersection'
                d = min(1,sum(histogramDistance(h(:,:,c,:,r+dr),h(:,:,c,:,r),'intersection'),3));
                d = 1-d;
            case 'chi2'
                d = min(1,sum(histogramDistance(h(:,:,c,:,r+dr),h(:,:,c,:,r),'chi2'),3));
                d = 1-d;
            case 'chi2-gaussian'
                d = min(1,sum(histogramDistance(h(:,:,c,:,r+dr),h(:,:,c,:,r),'chi2-gaussian',0.2),3));
                d = 1-d;
            case 'reconstruction'
                if channelIndex == 2
                    d = imageError(imgLab(:,:,c(1)),m(:,:,c(1),r),filters(r),et) + ...
                        imageError(imgLab(:,:,c(2)),m(:,:,c(2),r),filters(r),et);
                else
                    d = imageError(imgLab(:,:,c(1)),m(:,:,c(1),r),filters(r),et);
                end
            case 'sum'
                d = min(1, ...
                    sum(histogramDistance(h(:,:,1:3,:,r+dr),h(:,:,1:3,:,r),'chi2-gaussian',.2),3)...
                    + histogramDistance(h(:,:,4,:,r+dr),h(:,:,4,:,r),'chi2'));
            case 'combined'
                dmaxim = 1-min(1,sum(histogramDistance(h(:,:,1:3,:,r+dr),h(:,:,1:3,:,r),'chi2-gaussian',.2),3)...
                    + histogramDistance(h(:,:,4,:,r+dr),h(:,:,4,:,r),'chi2'));
                drecon = imageError(imgLab(:,:,1),m(:,:,1,r),filters(r),et) + ...
                    imageError(imgLab(:,:,2),m(:,:,2,r),filters(r),et) + ...
                    imageError(imgLab(:,:,3),m(:,:,3,r),filters(r),et);
                d = min(1,drecon + dmaxim);
            otherwise, error('Distance type is not supported')
        end
        % Disable annoying docking error that clutters the command line
        if strcmp(fh.WindowStyle, 'docked')
            warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
        end
        % Display image and then restore original patch
        subplot(122); imagesc(d,[0,1]); axis off image;
        title(sprintf('r_i=%d, r_o=%d, #bins=%d, %s, %s',...
                r,r+dr,B,channelType{channelIndex},distanceType{distanceIndex}));
        drawnow;
    end

    function changeRadius(fh,callbackData)
        r = r + callbackData.VerticalScrollCount; 
        if r < 1
            warning('Exceeded minimum radius (1)')
            r = 1;
        elseif r + dr > R
            warning('Exceeded maximum radius (%d). Decrease r or dr.',R)
            r = R-dr;
        else
            drawHistogramGradients(fh);
        end
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
            h = imageEncoding(imgBinned,filters,'hist',B);
        end
        drawHistogramGradients(fh)
    end

end