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
if ismatrix(imgRGB)
    channelType   = {'luminance','texture'};
else
    channelType   = {'luminance','color','texture'};
end
distanceType  = {'chi2','chi2-gaussian','reconstruction','sum','combined'};
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
mlab = imageEncoding(imgLab,filters,'average');
mrgb = imageEncoding(imgRGB,filters,'average');

% Plot figure and set callbacks
fh = figure; subplot(121); imshow(imgRGB); 
title(sprintf('Input image (%dx%d)',H,W));
set(fh, 'WindowButtonDownFcn',   @changeMethod);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
drawHistogramGradients(fh); % first draw

    function drawHistogramGradients(fh)
        % Compute distance
        et = 'mse';
%         dr = ceil(r/4);
        dr = ceil(r/(2+sqrt(6)));
        if strcmp(channelType{channelIndex},'color'), c = [2;3]; else c = channelIndex; end
        switch distanceType{distanceIndex}
            case {'intersection','chi2','chi2-gaussian'}
                h1 = h(:,:,c,:,r);
                h2 = h(:,:,c,:,r+dr)-h(:,:,c,:,r);
                h1 = bsxfun(@rdivide,h1,sum(h1,4));
                h2 = bsxfun(@rdivide,h2,sum(h2,4));
                d = sum(histogramDistance(h2,h1,distanceType{distanceIndex},0.2),3);
                d = 1-min(1,d);
            case 'reconstruction'
                if strcmp(et,'dssim')
                    d = imageError(imgRGB, mrgb(:,:,:,r),filters(r),'dssim');
                else
                    if ismatrix(imgRGB)
                        d = imageError(imgLab(:,:,1),mlab(:,:,1,r),filters(r),et);
                    else
                        d = (imageError(imgLab(:,:,1),mlab(:,:,1,r),filters(r),et) +...
                             imageError(imgLab(:,:,2),mlab(:,:,2,r),filters(r),et) +...
                             imageError(imgLab(:,:,3),mlab(:,:,3,r),filters(r),et))/2;
                    end
                    d = min(1,d); 
                end
            case 'sum'
                h1 = h(:,:,:,:,r);
                h2 = h(:,:,:,:,r+dr)-h(:,:,:,:,r);
                h1 = bsxfun(@rdivide,h1,sum(h1,4));
                h2 = bsxfun(@rdivide,h2,sum(h2,4));                
                d = sum(histogramDistance(h2(:,:,1:end-1,:),h1(:,:,1:end-1,:),'chi2-gaussian',0.2),3)...
                    + histogramDistance(h2(:,:,end,:),h1(:,:,end,:),'chi2');
                d = 1-min(1,d);
            case 'combined'
                h1 = h(:,:,:,:,r);
                h2 = h(:,:,:,:,r+dr)-h(:,:,:,:,r);
                h1 = bsxfun(@rdivide,h1,sum(h1,4));
                h2 = bsxfun(@rdivide,h2,sum(h2,4));                
                dmaxim = sum(histogramDistance(h2(:,:,1:end-1,:),h1(:,:,1:end-1,:),'chi2-gaussian',0.2),3)...
                    + histogramDistance(h2(:,:,end,:),h1(:,:,end,:),'chi2');
                dmaxim = 1-min(1,dmaxim);
                if strcmp(et,'dssim')
                    drecon = imageError(imgRGB,mrgb(:,:,:,r),filters(r),et);
                else
                    if ismatrix(imgRGB)
                        drecon = imageError(imgLab(:,:,1),mlab(:,:,1,r),filters(r),et);
                    else
                        drecon = (imageError(imgLab(:,:,1),mlab(:,:,1,r),filters(r),et) +...
                                  imageError(imgLab(:,:,2),mlab(:,:,2,r),filters(r),et) +...
                                  imageError(imgLab(:,:,3),mlab(:,:,3,r),filters(r),et))/2;
                    end
                end
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
            imgBinned = cat(3, binImage(imgLab,B),textonMap(imgRGB, B));
            h = imageEncoding(imgBinned,filters,'hist',B);
        end
        drawHistogramGradients(fh)
    end

end