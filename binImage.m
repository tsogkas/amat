function img = binImage(img,B)
% BINIMAGE Bins image from [0,1] into a set of discrete values [1,B].
img = max(1, ceil(bsxfun(@times, img, reshape(B,1,1,[]))));
