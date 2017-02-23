%% Google
google = imresize(imread('google.jpg'),0.5);
leopard= imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/134052.jpg'),0.5);
matgoogle = amat(google);
fprintf('Nonzero pixels: %d/%d (%.2f%%)\n', nnz(any(matgoogle.axis,3)), ...
    numel(matgoogle.radius),nnz(any(matgoogle.axis,3))/numel(matgoogle.radius)*100)

%% live demo
mat = amat(google,40,1e-3,1);

%% Texture image (leopard)
matleopard= amat(leopard);

%% Smoothed image
smoothed = L0Smoothing(leopard,2e-2,1.5);
imshow(smoothed)
matleopardSmoothed = amat(smoothed);

%% deer
deer = imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'),0.5,'bilinear');
matdeer = amat(deer);
deerSmoothed = L0Smoothing(deer);
matdeersmoothed = amat(deerSmoothed);

%% Boundary extraction
d = matgoogle.depth;
d = d/max(d(:));
imshow(1-d); title('Boundary strength')
imshow(bwmorph(matgoogle.depth<=8,'thin',inf))


