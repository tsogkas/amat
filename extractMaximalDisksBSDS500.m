%% Default directories and first image for testing data extraction
imageDir = '/home/tsogkas/datasets/BSDS500/images/train';
gtDir = '/home/tsogkas/datasets/BSDS500/groundtruth/train';
imageFiles = dir(fullfile(imageDir,'*.jpg'));
gtFiles = dir(fullfile(gtDir,'*.mat'));

% Plot obtained maximal circles as proof of concept
img = imread(fullfile(imageDir,imageFiles(1).name));
gt  = load(fullfile(gtDir, gtFiles(1).name)); gt = gt.groundTruth;

% Parameters
T = 50; % skeleton confidence threshold
minimumSegmentArea = 20;
resizeFactor = 0.5;
B = 32; % #bins
R = 40; % #scales
visualize = 0;

% Construct disk and ring filters
diskf = cell(R,1); for r=1:R, diskf{r} = disk(r); end
ringf = cell(R,1); for r=1:R, dr = ceil(r/10); ringf{r} = ring(r,r+dr); end


%% Compute binary MAT for all images
nImages = numel(imageFiles);
MAT(1:nImages) = struct('points',[],'radii',[],'img',[]);

% For all images
ticStart = tic;
for i=1:1
    img = imread(fullfile(imageDir,imageFiles(i).name));
    gt  = load(fullfile(gtDir, gtFiles(i).name)); gt = gt.groundTruth;
    img = imresize(img,resizeFactor, 'bilinear');
    [H,W,~] = size(img);
    nSegmentations = numel(gt);
    MAT(i).points = zeros(H,W,nSegmentations);
    MAT(i).radii  = zeros(H,W,nSegmentations);
    % For all segmentations
    for s=1:numel(gt)
        seg = imresize(gt{s}.Segmentation, resizeFactor, 'nearest');
        nSegments = numel(unique(seg));
        pmap = zeros(H,W);
        rmap = zeros(H,W);
        % For all segments in segmentation
        for j=1:nSegments
            segment = seg == j;
            if nnz(segment) >=  minimumSegmentArea
                [skel,r] = skeleton(segment);
                pmap = max(pmap,skel);
                rmap = max(rmap,sqrt(r)); % skeleton() returns squared distance for some reason
            end
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


%% Collect training examples
rng(0);
for i=1:nImages
    % First threshold the skeleton confidence to keep only medial points
    MAT(i).points = MAT(i).points > T;
    MAT(i).radii(~MAT(i).points) = 0;
    
    % Compute image features at all possible locations
    imrgb = MAT(i).img;
    imlab = rgb2labNormalized(imrgb);
    hdisk = imageEncoding(binImage(imlab,B),diskf, 'hist-normalized',B); 
    hring = imageEncoding(binImage(imlab,B),ringf, 'hist-normalized',B); 
    tdisk = imageEncoding(textonMap(imrgb,B),diskf,'hist-normalized',B);
    tring = imageEncoding(textonMap(imrgb,B),ringf,'hist-normalized',B);
    
    % Get positive samples (true maximal disks). Not-thinning the skeleton
    % maps allows for some slack in the positive samples considered.
    
    % Get negative samples, centered at the same points as positive
    % samples, but at different scales, and create ordering.
    
    % Get random negative samples from other points in the image.
    
end




