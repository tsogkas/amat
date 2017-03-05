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
% matrefined = refineMAT(mat);
% The group labels are already sorted and first label is zero (background)
%%
numBranches = max(mat.branches(:)); 
branchOld   = bsxfun(@eq, mat.branches, reshape(1:numBranches,1,1,[]));
% sort branches
dilated     = false(size(branchOld));
isodilated  = false(size(branchOld));
thindil     = false(size(branchOld));
thinisodil  = false(size(branchOld));
skeldil     = false(size(branchOld));
skelisodil  = false(size(branchOld));
cover       = false(size(branchOld));
thincover   = false(size(branchOld));
skelcover   = false(size(branchOld));
for i=1:numBranches
    cover(:,:,i) = mat2mask(mat.radius.*branchOld(:,:,i), mat.scales)>0;
end
[~,idxSorted] = sort(sum(sum(cover)),'descend');
cover = cover(:,:,idxSorted);
branchOld = branchOld(:,:,idxSorted);
r = 3;
SE = strel('disk',r);
for i=1:numBranches
    dilated(:,:,i)    = imdilate(branchOld(:,:,i),SE);
    isodilated(:,:,i) = bwdist(branchOld(:,:,i)) <= r;
    thindil(:,:,i)    = bwthin(dilated(:,:,i));
    thinisodil(:,:,i) = bwthin(isodilated(:,:,i));
    thincover(:,:,i)  = bwthin(cover(:,:,i));
    skelcover(:,:,i)  = bwmorph(cover(:,:,i),'skel',inf);
end

%%
for i=1:min(numBranches,10)
    figure, imshowpair(cover(:,:,i),branchOld(:,:,i)); 
    title('Branch and cover')
    figure, imshowpair(thindil(:,:,i),thinisodil(:,:,i)); 
    title(sprintf('Thinned dilation and isotropic dilation (r=%d)',r))
    figure, imshow(skelcover(:,:,i)); 
    title('Skeletonized cover')
%     figure, imshowpair(thincover(:,:,i),skelcover(:,:,i)); 
%     title('Thinned and skeletonized cover')
end

%% Edge boxes (PACSCAL images)
pascalpath = '/home/tsogkas/datasets/VOC2007/VOCdevkit/VOC2007/JPEGImages/';
imgfile = '000025';
img = imresize(imread(fullfile(pascalpath,[imgfile, '.jpg'])),0.5,'bilinear');
smoothed = L0Smoothing(img);
mat = amat(smoothed);
mat.branches = groupMedialPoints(mat);
%%
[seg,segments] = mat2seg(mat,0.9);
e = seg2edges(seg);
