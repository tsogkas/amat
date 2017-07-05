function out = imcropCenter(A, height, width)
%IMCROPCENTER crop image A to size height*width

[H,W,~] = size(A);
if (height > H || width > W)
    error('height and width must be less than the size of A.')
end
pad = [H - height, W - width] / 2;
out = A(ceil(pad(1))+1:H-floor(pad(1)), ceil(pad(2))+1:W-floor(pad(2)), :);

end

