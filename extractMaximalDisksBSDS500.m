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
slack = 3; % 
visualize = 0;

% Construct disk and ring filters
diskf = cell(R,1); for r=1:R, diskf{r} = disk(r); end
ringf = cell(R,1); for r=1:R, dr = ceil(r/10); ringf{r} = ring(r,r+dr); end


%% Compute binary MAT for all images
nImages = numel(imageFiles);
MAT(1:nImages) = struct('pts',[],'rad',[],'img',[],'seg',[],'bnd',[]);

% For all images
ticStart = tic;
for i=1:1
    img = imread(fullfile(imageDir,imageFiles(i).name));
    gt  = load(fullfile(gtDir, gtFiles(i).name)); gt = gt.groundTruth;
    img = imresize(img,resizeFactor, 'bilinear');
    [H,W,~] = size(img);
    nSegmentations = numel(gt);
    MAT(i).pts = zeros(H,W,nSegmentations);
    MAT(i).rad = zeros(H,W,nSegmentations);
    MAT(i).bnd = false(H,W,nSegmentations);
    % For all segmentations
    for s=1:nSegmentations
        seg = imresize(gt{s}.Segmentation, resizeFactor, 'nearest');
        nSegments = numel(unique(seg));
        pmap = zeros(H,W);
        rmap = zeros(H,W);
        bmap = false(H,W);
        % For all segments in segmentation
        for j=1:nSegments
            segment = seg == j;
            if nnz(segment) >=  minimumSegmentArea
                [skel,r] = skeleton(segment);
                bmap = bmap | bwperim(segment);
                pmap = max(pmap,skel);
                rmap = max(rmap,sqrt(r)); % skeleton() returns squared distance for some reason
            end
        end
        MAT(i).pts(:,:,s) = pmap;
        MAT(i).rad(:,:,s) = rmap;
        MAT(i).bnd(:,:,s) = bmap;
        MAT(i).seg(:,:,s) = seg;
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
rng(0); % reset rng for reproducibility 
for i=1:1
    % First threshold the skeleton confidence to keep only medial points
    p = MAT(i).pts;
    r = MAT(i).rad;
    p = p > T;
    r(~p) = 0;
    
    % Compute image features at all possible locations
    imrgb = MAT(i).img;
    imlab = rgb2labNormalized(imrgb);
    nbpts = imageEncoding(double(MAT(i).bnd), diskf, 'sum'); 
    hdisk = imageEncoding(binImage(imlab,B),diskf, 'hist-normalized',B); 
    hring = imageEncoding(binImage(imlab,B),ringf, 'hist-normalized',B); 
    tdisk = imageEncoding(textonMap(imrgb,B),diskf,'hist-normalized',B);
    tring = imageEncoding(textonMap(imrgb,B),ringf,'hist-normalized',B);
    
    % Create flag array of points that are used as training samples
    [H,W,~] = size(imrgb);
    used = zeros(H,W,R); 
    
    % Get positive samples (true maximal disks). Not-thinning the skeleton
    % maps allows for some slack in the positive samples considered.
    ind = find(p);
    [ypos,xpos,~] = ind2sub(size(p),ind);
    rpos = round(r(ind));
    rinvalid = rpos > R | rpos <= 0;
    ypos(rinvalid) = []; % throw away samples with radius > R or radius <=0
    xpos(rinvalid) = [];
    rpos(rinvalid) = [];
    pos = [ypos,xpos,rpos];
    pos = unique(pos,'rows'); % throw away duplicate samples
    ind = sub2ind([H,W,R], ypos,xpos,rpos);
    used(ind) = 1; % +1 for positive samples
    
    % LEAVE THAT FOR LATER 
    % Get "hard" negative samples, that have very small distance from
    % boundaries or touch boundaries only at one point. These can be
    % centered at the same points as positives but have different scales,
    % and should have a lower cost than disks that cross boundaries, or
    % that are fully contained in the interior of the shape.
       
    % Get random negative samples from other points in the image. Use a
    % balanced set of #neg == #pos.
    % TODO: Maybe change the distribution used to draw the radii for the
    % negative samples from uniform to one that favors smaller and medium
    % scales.
    nPos = size(pos,1);
    ind = find(used == 0);
    ind = ind(randsample(1:length(ind), nPos));
    [yneg,xneg,rneg] = ind2sub([H,W,R],ind);
    rinvalid = rneg > R | rneg <= 0;
    yneg(rinvalid) = []; % throw away samples with radius > R or radius <=0
    xneg(rinvalid) = [];
    rneg(rinvalid) = [];
    neg = [yneg,xneg,rneg];
    neg = unique(neg,'rows'); % throw away duplicate samples
    ind = sub2ind([H,W,R], yneg,xneg,rneg);
    used(ind) = -1; % -1 for negative samples
    
    if visualize
        figure(1);
        imshow(imrgb); step = 20;
        viscircles([xpos(1:step:end),ypos(1:step:end)],rpos(1:step:end),...
            'Color','g','EnhanceVisibility',true,'Linewidth',0.5);
        viscircles([xneg(1:step:end),yneg(1:step:end)],rneg(1:step:end),...
            'Color','r','EnhanceVisibility',true,'Linewidth',0.5);
    end
    
    % Build feature vectors for positive and negative samples
    % Concatenate color and texture histograms.
    hdisk = cat(3,hdisk,tdisk); % HxWxCxBxR (C = 4)
    hring = cat(3,hring,tring); % HxWxCxBxR (C = 4)
    [H,W,C,B,R] = size(hdisk);
    hdisk = reshape(permute(hdisk,[1,2,5,3,4]), H*W*R, []);
    hring = reshape(permute(hring,[1,2,5,3,4]), H*W*R, []);
    f = [hdisk; hring];
    fpos = f(used == 1,:);
    fneg = f(used == -1,:);
    
    
    
end




