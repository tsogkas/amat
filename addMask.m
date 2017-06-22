function out = addMask(I, y, x, msk)
% Add msk to image I centered at y,x

[H,W,~] = size(I);

% make sure fil has a center
msk = padarray(msk, double(mod(size(msk), 2) == 0), 'pre');

[ry,rx,~] = size(msk);
ry = floor(ry / 2);
rx = floor(rx / 2);

% Since we are using convolution, the image should be rotated by 180 deg.
msk = imrotate(msk,180);

y1 = y - ry; y2 = y + ry; x1 = x - rx; x2 = x + rx;

if (y-ry < 1)
    msk = msk(2+ry-y:end,:,:);
    y1 = 1;
end
if (x-rx < 1)
    msk = msk(:,2+rx-x:end,:);
    x1 = 1;
end
if (y+ry > H)
    msk = msk(1:end-(y+ry-H),:,:);
    y2 = H;
end
if (x+rx > W)
    msk = msk(:,1:end-(x+rx-W),:);
    x2 = W;
end

I(y1:y2,x1:x2,:) = I(y1:y2,x1:x2,:) + msk;
out = I;

end
