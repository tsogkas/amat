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
mat = matbird;
mat.branches = groupMedialPoints(mat);
branches = mat.branches;
radius   = mat.radius;
% Compute the depth contribution of each branch separately.
[H,W] = size(mat.depth);
numBranches = max(branches(:));
depthBranch = zeros(H,W,numBranches);
for i=1:numBranches
    depthBranch(:,:,i) = mat2mask(radius .* double(branches == i));
end

%%
% The edges are the points in the image that are covered by few disks
% belonging in the current branch AND by few disks belonging to OTHER
% branches. We use the edge strength of all the other branches to weigh the
% strength of the edges of the current branch.
d  = depthBranch;
dc = bsxfun(@minus,mat.depth,d); % complement of the depth for a branch
dsoft  = 1-bsxfun(@rdivide,d, max(eps,max(max(d))));
dcsoft = 1-bsxfun(@rdivide,dc,max(eps,max(max(dc))));
dpsoft = dsoft .* dcsoft;

branchImportance = sum(sum(depthBranch>0))/(size(depthBranch,1)*size(depthBranch,2));
[importanceSorted,idxSorted] = sort(branchImportance,'descend');
d = d(:,:,idxSorted);
dc= dc(:,:,idxSorted);
dsoft  = dsoft(:,:,idxSorted);
dcsoft = dcsoft(:,:,idxSorted);
dpsoft = dpsoft(:,:,idxSorted);
% e = prod(dp(:,:,idxSorted(1:10)));


%% Test refineMAT()
mat = mattotem;
mat.branches = groupMedialPoints(mat);
matrefined = refineMAT(mat);

