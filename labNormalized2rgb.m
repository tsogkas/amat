function imRGB = labNormalized2rgb(imLAB)
% LABNORMALIZED2RGB Convert input image from the normalized in [0,1] 
%   CIE-Lab to the RGB color space ([0,1]).
% 
%   imRGB = LABNORMALIZED2RGB(imLAB)
% 
% Stavros Tsogkas, <tsogkas@cs.toronto.edu>
% Last update: February 2017

imRGB = imLAB;
if ~ismatrix(imRGB) % RGB image or image stack
    % We use the min and max values for the Lab channels to scale in [0,1].
    % NOTE: The Berkeley pb code used the values abmin = -73 and abmax = 95
    abmin = -128; abmax = 127;
    imRGB(:,:,1) = imRGB(:,:,1) * 100;
    imRGB(:,:,2:3) = (abmax-abmin) * imRGB(:,:,2:3) + abmin;
    
    if verLessThan('matlab','R2014b')
        imRGB = applycform(imRGB, makecform('lab2srgb'));
    else
        imRGB = lab2rgb(imRGB);
    end
end