%% Google
google = imresize(imread('google.jpg'),0.5);
matgoogle = amat(google,40,1e-4);
fprintf('Nonzero pixels: %d/%d (%.2f%%)\n', nnz(any(matgoogle.axis,3)), ...
    numel(matgoogle.radius),nnz(any(matgoogle.axis,3))/numel(matgoogle.radius)*100)

%% live demo
% mat = amat(google,40,1e-3,1);

%% Texture image (leopard)
matleopard= amat(leopard);

%% Smoothed image
leopard= imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/134052.jpg'),0.5);
smoothed = L0Smoothing(leopard,2e-2,1.5);
imshow(smoothed)
matleopardSmoothed = amat(smoothed);

%% deer
deer = imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'),0.5,'bilinear');
% matdeer = amat(deer);
deerSmoothed = L0Smoothing(deer);
matdeer = amat(deerSmoothed);

%% totems
totem = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/101085.jpg'),0.5,'bilinear');
totemSmoothed = L0Smoothing(totem);
mattotem = amat(totemSmoothed,40,1e-4);

%% bird
bird = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/42049.jpg'),0.5,'bilinear');
birdSmoothed = L0Smoothing(bird);
matbird = amat(birdSmoothed,40,1e-4);

%% plane
plane = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/37073.jpg'),0.5,'bilinear');
planeSmoothed = L0Smoothing(plane);
matplane = amat(planeSmoothed);

%% seeds
seeds = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/58060.jpg'),0.5,'bilinear');
seedsSmoothed = L0Smoothing(seeds,3e-2);
matseeds = amat(seedsSmoothed);

%% vase
vase = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/227092.jpg'),0.5,'bilinear');
vaseSmoothed = L0Smoothing(vase);
matvase = amat(vaseSmoothed);

%% penguin
penguin = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/106024.jpg'),0.5,'bilinear');
penguinSmoothed = L0Smoothing(penguin);
matpenguin = amat(penguinSmoothed,40,1e-4);

%% manstick
manstick = imresize(imread('/home/tsogkas/datasets/BSDS500/images/val/101087.jpg'),0.5,'bilinear');
% manstick = imread('/home/tsogkas/datasets/BSDS500/images/val/101087.jpg');
manstickSmoothed = L0Smoothing(manstick);
matmanstick = amat(manstickSmoothed,40,1e-4);

%% Boundary extraction
mat = matmanstick;
mat.branches = groupMedialPoints(mat);

%%
minCoverage = 0.95;
seg = mat2seg(mat,minCoverage);
e = seg2edges(seg);
figure; imshow(label2rgb(seg));
figure; imshow(label2rgb(e));

%% Test refineMAT()
mat = mattotem;
mat.branches = groupMedialPoints(mat);
matrefined = refineMAT(mat);

