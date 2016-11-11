function imLAB = rgb2lab(imRGB)
% RGB2LAB Convert input image from RGB to the CIE-Lab color space
% 
%   imLAB = RGB2LAB(imRGB)
% 
% Stavros Tsogkas, <tsogkas@cs.toronto.edu>
% Last update: October 2016 

if ismatrix(imRGB) % grayscale image
    imLAB = imRGB;
else % RGB image
    % convert gamma-corrected image to LAB and scale values into [0,1]
    % min and max values for a,b channels of LAB
    % used to scale values into the unit interval
    abmin = -73;
    abmax = 95;
    % CHECK THIS: DO WE NEED GAMMA CORRECTION WITH makecform??
%   gamma = 2.5; lab = RGB2Lab(im.^gamma); 
    imLAB = applycform(imRGB, makecform('srgb2lab'));
    imLAB(:,:,1) = imLAB(:,:,1) ./ 100;
    imLAB(:,:,2) = (imLAB(:,:,2) - abmin) ./ (abmax-abmin);
    imLAB(:,:,3) = (imLAB(:,:,3) - abmin) ./ (abmax-abmin);
    imLAB(:,:,2) = max(0,min(1,imLAB(:,:,2)));
    imLAB(:,:,3) = max(0,min(1,imLAB(:,:,3)));
end