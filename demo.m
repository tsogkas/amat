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
% pascalpath = '/home/tsogkas/datasets/VOC2007/VOCdevkit/VOC2007/JPEGImages/';
% imgfile = '000025';
% img = imresize(imread(fullfile(pascalpath,[imgfile, '.jpg'])),0.5,'bilinear');

% Load model and default edge Boxes options
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;
opts = edgeBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% Read image and compute mat, branches, and edges
img = imread('peppers.png');
smoothed = L0Smoothing(imresize(img,0.5,'bilinear'));
mat = amat(smoothed);
mat.branches = groupMedialPoints(mat);
[seg,segments] = mat2seg(mat,0.9);
eg = single(seg2edges(seg)); % labelled edges
es = single(bwthin(imresize(eg>0, [size(img,1), size(img,2)],'nearest'))); % binarized edges
O2  = edgeOrient(es);
[E,O]=edgesDetect(img,model); E=edgesNmsMex(E,O,2,0,1,model.opts.nThreads);
%% First test if it's working with our binarized edge output
o = opts;
bbs1 = edgeBoxes(img,model,opts);
bbs2 = edgeBoxesMex(single(es),O2,o.alpha,o.beta,o.eta,o.minScore,o.maxBoxes,...
  o.edgeMinMag,o.edgeMergeThr,o.clusterMinMag,...
  o.maxAspectRatio,o.minBoxArea,o.gamma,o.kappa);



%%
gt=[122 248 92 65; 193 82 71 53; 410 237 101 81; 204 160 114 95; ...
  9 185 86 90; 389 93 120 117; 253 103 107 57; 81 140 91 63];
if(0), gt='Please select an object box.'; disp(gt); figure(1); imshow(I);
  title(gt); [~,gt]=imRectRot('rotate',0); gt=gt.getPos(); end
gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs1),.7);
figure(1); bbGt('showRes',img,gtRes,dtRes(dtRes(:,6)==1,:));
title('green=matched gt  red=missed gt  dashed-green=matched detect');

[gtRes,dtRes]=bbGt('evalRes',gt,double(bbs2),.7);
figure(2); bbGt('showRes',img,gtRes,dtRes(dtRes(:,6)==1,:));
title('green=matched gt  red=missed gt  dashed-green=matched detect');
