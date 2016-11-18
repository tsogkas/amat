function imLAB = rgb2labNormalized(imRGB)
% RGB2LABNORMALIZED Convert input image from RGB to the CIE-Lab color space
%   and normalize in the unit interval [0,1].
% 
%   imLAB = RGB2LABNORMALIZED(imRGB)
% 
% Stavros Tsogkas, <tsogkas@cs.toronto.edu>
% Last update: November 2016 

if ismatrix(imRGB) % grayscale image
    imLAB = imRGB;
else % RGB image
    if verLessThan('matlab','R2014b')
        imLAB = applycform(imRGB, makecform('srgb2lab'));
    else
        imLAB = rgb2lab(imRGB);
    end
    % We use the min and max values for the Lab channels to scale in [0,1].
    % NOTE: The Berkeley pb code used the values abmin = -73 and abmax = 95
    abmin = -128; abmax = 127;
    imLAB(:,:,1) = imLAB(:,:,1) ./ 100;
    imLAB(:,:,2) = (imLAB(:,:,2) - abmin) ./ (abmax-abmin);
    imLAB(:,:,3) = (imLAB(:,:,3) - abmin) ./ (abmax-abmin);
    
    % Make sure we are in [0,1]
    imLAB = max(0, min(1, imLAB));
end