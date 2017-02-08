function BMAX500 = constructBMAX500(varargin)
% Build medial axis transform dataset based on the BSDS500

opts = {'resize',             1,... % this can be a HxW vector
        'skelThresh',         50,...% high --> aggressive skeleton pruning
        'minRadius',          3,... % we ignore points with radius < minRadius
        'parpoolSize',        inf   % set to 0 to run serially 
       };
opts = parseVarargin(opts, varargin, 'struct');


% Compute skeletons for all available segmentations and segments in BSDS500
paths = setPaths();
savePath = fullfile(paths.amat.output, ['BMAX500-resize' num2str(opts.resize) '.mat']);

try
    disp('Loading BMAX500 from disk...')
    BMAX500 = load(savePath); BMAX500 = BMAX500.BMAX500;
catch
    disp('File was not found. Computing medial axes for BSDS500 dataset...')
    BMAX500 = struct('train',[],'val',[],'test',[],'opts',opts);
    % For all subsets of BSDS500
    for set = {'train','val','test'}
        % Setup image and groundtruth paths
        imDir = fullfile(paths.bsds500im, set{1});
        gtDir = fullfile(paths.bsds500gt, set{1});
        imFiles = dir(fullfile(imDir,'*.jpg'));
        gtFiles = dir(fullfile(gtDir,'*.mat'));
        assert(numel(imFiles) == numel(gtFiles), '#img ~= #seg')
        nImages = numel(imFiles);
        matgt = repmat(struct(...
            'pts',[],'rad',[],'img',[],'seg',[],'bnd',[],'iid',[]), [1,nImages]);
        % For all images
        ticStart = tic;
        parfor (i=1:nImages, opts.parpoolSize)
            img = imread(fullfile(imDir,imFiles(i).name));
            gt  = load(fullfile(gtDir, gtFiles(i).name)); gt = gt.groundTruth;
            if numel(opts.resize) == 2 || opts.resize ~= 1
                img = imresize(img,opts.resizeFactor, 'bilinear');
            end
            [H,W,~] = size(img);
            nSegmentations = numel(gt);
            [~,iid] = fileparts(imFiles(i).name);
            matgt(i).iid = iid;   
            matgt(i).img = img;   
            matgt(i).pts = zeros(H,W,nSegmentations,'uint8');
            matgt(i).rad = zeros(H,W,nSegmentations,'uint8');
            matgt(i).seg = zeros(H,W,nSegmentations,'uint8');
            matgt(i).bnd = false(H,W,nSegmentations);
            % For all segmentations
            for s=1:nSegmentations
                seg = gt{s}.Segmentation;
                if numel(opts.resize) == 2 || opts.resize ~= 1
                    seg = imresize(seg, opts.resize, 'nearest');
                end
                nSegments = numel(unique(seg));
                pmap = zeros(H,W);
                rmap = zeros(H,W);
                bmap = false(H,W);
                % Compute skeletons of all segments
                for j=1:nSegments
                    segment = seg == j;
                    [skel,r] = skeleton(segment);
                    bmap = bmap | bwperim(segment);
                    pmap = max(pmap,skel);
                    rmap = max(rmap,sqrt(r)); % skeleton() returns squared distance for some reason
                end
                matgt(i).pts(:,:,s) = pmap;
                matgt(i).rad(:,:,s) = rmap;
                matgt(i).bnd(:,:,s) = bmap;
                matgt(i).seg(:,:,s) = seg;             
            end
            progress(['Computing MAT (' set{1} ')...'],i,nImages,ticStart,-1);
        end
        BMAX500.(set{1}) = matgt;
    end
    save(savePath, 'BMAX500');
end

% Refine dataset 
disp('Refining dataset')
for set = {'train','val','test'}
    if ~isfield(BMAX500, set{1}), continue; end
    matgt = BMAX500.(set{1});
    for i=1:numel(BMAX500.(set{1}))
        % Prune points with low confidence
        matgt(i).pts = matgt(i).pts >= opts.skelThresh;
        % Prune points with very small radius
        matgt(i).pts(matgt(i).rad < opts.minRadius) = false;
    end
    BMAX500.(set{1}) = matgt;
end




