%% Default directories and first image for testing data extraction

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
paths = setPaths();
imDir = paths.bsds500imTrain; imFiles = dir(fullfile(imDir,'*.jpg'));
gtDir = paths.bsds500gtTrain; gtFiles = dir(fullfile(gtDir,'*.mat'));
assert(numel(imFiles) == numel(gtFiles), '#img ~= #seg')

savePath = fullfile(paths.output, 'MATBSDS500TRAIN.mat');
% nImages = numel(imFiles);
nImages = 1;
MAT(1:nImages) = struct('pts',[],'rad',[],'img',[],'seg',[],'bnd',[]);

try
    disp('Loading extracted MAT for BSDS500 training images...')
    MAT = load(savePath); MAT = MAT.MAT;
catch
    disp('File was not found. Computing MAT for BSDS500 train images...')
    % For all images
    ticStart = tic;
    for i=1:nImages
        img = imread(fullfile(imDir,imFiles(i).name));
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
        progress('Extracting MAT for BSDS500 training images...',i,nImages,ticStart,-1);
    end
    save(savePath, 'MAT');
end

%% Collect training examples

savePath = fullfile(paths.output, 'FEATBSDS500TRAIN.mat');
try
    disp('Loading features extracted from BSDS500 training images...')
    tmp = load(savePath); fposall = tmp.fposall; fnegall = tmp.fnegall; clear tmp;
catch
    disp('File was not found. Computing features for BSDS500 training images...')
    rng(0); % reset rng for reproducibility 
    fposall = []; % features for positives from all images will be stored here
    fnegall = []; % features for negatives from all images will be stored here
    ticStart = tic;
    for i=1:nImages
        % First threshold the skeleton confidence to keep only medial points
        p = MAT(i).pts;
        r = MAT(i).rad;
        p = p > T;
        r(~p) = 0;

        % Compute image features at all possible locations
        imrgb = MAT(i).img;
        imlab = rgb2labNormalized(imrgb);
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
        pos = unique([xpos,ypos,rpos],'rows'); % throw away duplicate samples
        xpos = pos(:,1); ypos = pos(:,2); rpos = pos(:,3);
        ind = sub2ind([H,W,R], ypos,xpos,rpos);
        used(ind) = 1; % +1 for positive samples

        % LEAVE THAT FOR LATER 
        % Get "hard" negative samples, that have very small distance from
        % boundaries or touch boundaries only at one point. These can be
        % centered at the same points as positives but have different scales,
        % and should have a lower cost than disks that cross boundaries, or
        % that are fully contained in the interior of the shape.

        % Get random negative samples from other points in the image. Use a
        % balanced set of #neg Some tools below are only available to our subscribers or users with an online account.== #pos.
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
        neg = unique([xneg,yneg,rneg],'rows'); % throw away duplicate samples
        xneg = neg(:,1); yneg = neg(:,2); rneg = neg(:,3);
        ind = sub2ind([H,W,R], yneg,xneg,rneg);
        used(ind) = -1; % -1 for negative samples

        if visualize
            figure(1);
            imshow(imrgb); step = 20;
            viscircles([xpos(1:step:end),ypos(1:step:end)],rpos(1:step:end),'Color','g','Linewidth',0.5);
            viscircles([xneg(1:step:end),yneg(1:step:end)],rneg(1:step:end),'Color','r','Linewidth',0.5);
        end

        % Build feature vectors for positive and negative samples
        % Merge color and texture histograms.
        hdisk = cat(3,hdisk,tdisk); % HxWxCxBxR (C = 4)
        hring = cat(3,hring,tring); % HxWxCxBxR (C = 4)
        hdisk = reshape(permute(hdisk,[1,2,5,3,4]), H*W*R, []);
        hring = reshape(permute(hring,[1,2,5,3,4]), H*W*R, []);
        clear tdisk tring

        % Concatenate all features: 
        % 1) Histograms corresponding to disks
        % 2) HisSome tools below are only available to our subscribers or users with an online account.tograms corresponding to the outer rings
        % 3) Chi-squared distance terms (h1-h2).^2 ./ (h1+h2+eps)
        % 4) Fixed term 1/r used as "scale feature"
        fr = repmat(reshape((1:R)./R,1,1,[]), [H,W,1]); % scale feature    
%         fall = [hdisk, hring, (hdisk-hring).^2 ./ (hdisk+hring+eps), fr(:)];
        fall = [hdisk, hring, (hdisk-hring).^2];
        clear hdisk hring fr;
        fpos = fall(used == +1,:); % nPos x 2*B 
        fneg = fall(used == -1,:); 
        clear fall;

        % Add to the feature matrix
        fposall = [fposall; fpos];
        fnegall = [fnegall; fneg];
        progress('Extracting features for BSDS500 training images...',i,nImages,ticStart);
    end
    save(savePath,'-v7.3','fposall','fnegall');
end

%% Standard SVM

% TODO: Should I add the bias term myself?
global X
X = [fposall; fnegall];
Y = [ones(size(fposall,1),1); -ones(size(fnegall,1),1)];

% Set training options and train RankSVM
linear = 1;
lambda = 1e-3;
[w,b,obj] = primal_svm(linear,Y,lambda);

scores = X*w + b; 
% nnz(scores(Y==1)>0) / nnz(Y==1)
% nnz(scores(Y==-1)>0)/nnz(Y==-1)
nnz(scores(Y==-1)>min(scores(Y==1)))



%% Train RankSVM

% Create pairs that are used as training examples in the RankSVM
% TODO: Should I add the bias term myself?
global X
X = fposall - fnegall;
fnegall = fnegall(randperm(size(fnegall,1)),:);
Y = ones(size(X,1),1);
X = [X; fnegall - fposall];
Y = [Y; -Y];

% Set training options and train RankSVM
linear = 1;
lambda = 1e-9;
[w,b,obj] = primal_svm(linear,Y,lambda);

pscores = fposall*w + b; 
nscores = fnegall*w + b; 

nnz(nscores>min(pscores))



