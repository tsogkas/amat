function imshow3(im1,im2,im3)
% IMSHOW3 Displays three images side by side, using the IMSHOW() function.
% 
%   IMSHOW3(im1,im2,im3) acts as a shortcut of IMSHOW ans SUBPLOT combinations
%   to display images im1 and im2 side by side.
% 
%   IMSHOW3({im1,title1},{im2,title2},{im3,title3}) adds an optional title for each
%   image
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: October 2016

title1 = [];
title2 = [];
title3 = [];
if iscell(im1)
    if numel(im1) == 2, title1 = im1{2}; end
    im1 = im1{1};
end
if iscell(im2)
    if numel(im2) == 2, title2 = im2{2}; end
    im2 = im2{1};
end
if iscell(im3)
    if numel(im3) == 2, title3 = im3{2}; end
    im3 = im3{1};
end

subplot(131); imshow(im1); title(title1);
subplot(132); imshow(im2); title(title2);
subplot(133); imshow(im3); title(title3);