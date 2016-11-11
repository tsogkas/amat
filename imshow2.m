function imshow2(im1,im2)
% IMSHOW2 Displays two images side by side, using the IMSHOW() function.
% 
%   IMSHOW2(im1,im2) acts as a shortcut of IMSHOW ans SUBPLOT combinations
%   to display images im1 and im2 side by side.
% 
%   IMSHOW2({im1,title1},{im2,title2}) adds an optional title for each
%   image
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: October 2016

title1 = [];
title2 = [];
if iscell(im1)
    if numel(im1) == 2, title1 = im1{2}; end
    im1 = im1{1};
end
if iscell(im2)
    if numel(im2) == 2, title2 = im2{2}; end
    im2 = im2{1};
end

subplot(121); imshow(im1); title(title1);
subplot(122); imshow(im2); title(title2);