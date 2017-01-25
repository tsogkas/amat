%% Default directories and first image for testing data extraction
imageDir = '/home/tsogkas/datasets/BSDS500/images/train';
gtDir = '/home/tsogkas/datasets/BSDS500/groundtruth/train';
imageFiles = dir(fullfile(imageDir,'*.jpg'));
gtFiles = dir(fullfile(gtDir,'*.mat'));

% Plot obtained maximal circles as proof of concept
img = imread(fullfile(imageDir,imageFiles(1).name));
gt  = load(fullfile(gtDir, gtFiles(1).name)); gt = gt.groundTruth;
[H,W,C] = size(img);
T = 35; % use a higher threshold to make sure we get good disks
visualize = 0;

%% Compute binary MAT for all images
nImages = numel(imageFiles);
MAT(1:nImages) = struct('points',[],'radii',[],'img',[]);

% For all images
ticStart = tic;
parfor i=1:nImages
    img = imread(fullfile(imageDir,imageFiles(i).name));
    gt  = load(fullfile(gtDir, gtFiles(i).name)); gt = gt.groundTruth;
    nSegmentations = numel(gt);
    MAT(i).points = zeros(H,W,nSegmentations);
    MAT(i).radii  = zeros(H,W,nSegmentations);
    % For all segmentations
    for s=1:numel(gt)
        seg = gt{s}.Segmentation;
        nSegments = numel(unique(seg));
        pmap = zeros(H,W);
        rmap = zeros(H,W);
        % For all segments in segmentation
        for j=1:nSegments
            segment = seg == j;
            [skel,r] = skeleton(segment);
            p = bwmorph(skel > T,'skel','inf');
            pmap(p) = skel(p);
            rmap(p) = sqrt(r(p)); % skeleton() returns the squared distance for some reason
        end
        MAT(i).points(:,:,s) = pmap;
        MAT(i).radii(:,:,s)  = rmap;
        MAT(i).img = img;
        
        % Visualize
        if visualize
            figure(s); imagesc(seg); axis image off;
            [yy,xx] = find(pmap > 0);
            rr = rmap(pmap > 0);
            assert(numel(xx) == numel(yy) && numel(xx) == numel(rr))
            viscircles([xx,yy],rr,'Color','r','EnhanceVisibility',true,'Linewidth',0.5);
        end
    end
    progress('Extracting MAT for training images...',i,nImages,ticStart,-1);
end

%% Extract training examples
