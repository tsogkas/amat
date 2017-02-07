function BMAX500 = constructBMAX500(varargin)
% Build medial axis transform dataset based on the BSDS500

opts = {'minSegmentPixels',   50,...
        'resizeFactor',       1,...
        'visualize',          false
       };
opts = parseVarargin(opts, varargin, 'struct');   

paths = setPaths();
savePath = fullfile(paths.output, 'BMAX500.mat');
try
    disp('Loading BMAX500 from disk...')
    BMAX500 = load(savePath); BMAX500 = BMAX500.BMAX500;
catch
    disp('File was not found. Computing medial axes for BSDS500 dataset...')
    BMAX500 = struct('train',[],'val',[],'test',[],'opts',opts);
    % For all subsets of BSDS500
    for set = {'train','val','test'}
        disp(['Extracting BMAX500, ''' set{1} ''' set...'])
        % Setup image and groundtruth paths
        imDir = fullfile(paths.bsds500im, set{1}); 
        gtDir = fullfile(paths.bsds500gt, set{1}); 
        imFiles = dir(fullfile(imDir,'*.jpg'));
        gtFiles = dir(fullfile(gtDir,'*.mat'));
        assert(numel(imFiles) == numel(gtFiles), '#img ~= #seg')
        nImages = numel(imFiles);
        matgt(1:nImages) = struct('pts',[],'rad',[],'img',[],'seg',[],'bnd',[]);
        % For all images
        ticStart = tic;
        for i=1:nImages
            img = imread(fullfile(imDir,imFiles(i).name));
            gt  = load(fullfile(gtDir, gtFiles(i).name)); gt = gt.groundTruth;
            if numel(opts.resizeFactor) == 2 || opts.resizeFactor ~= 1
                img = imresize(img,opts.resizeFactor, 'bilinear');
            end
            [H,W,~] = size(img);
            nSegmentations = numel(gt);
            matgt(i).pts = zeros(H,W,nSegmentations);
            matgt(i).rad = zeros(H,W,nSegmentations);
            matgt(i).bnd = false(H,W,nSegmentations);
            % For all segmentations
            for s=1:nSegmentations
                seg = gt{s}.Segmentation;
                if numel(opts.resizeFactor) == 2 || opts.resizeFactor ~= 1
                    seg = imresize(seg, opts.resizeFactor, 'nearest');
                end
                nSegments = numel(unique(seg));
                pmap = zeros(H,W);
                rmap = zeros(H,W);
                bmap = false(H,W);
                % For all segments in segmentation
                for j=1:nSegments
                    segment = seg == j;
                    if nnz(segment) >= opts.minSegmentPixels
                        [skel,r] = skeleton(segment);
                        bmap = bmap | bwperim(segment);
                        pmap = max(pmap,skel);
                        rmap = max(rmap,sqrt(r)); % skeleton() returns squared distance for some reason
                    end
                end
                matgt(i).pts(:,:,s) = pmap;
                matgt(i).rad(:,:,s) = rmap;
                matgt(i).bnd(:,:,s) = bmap;
                matgt(i).seg(:,:,s) = seg;
                matgt(i).img = img;
                
                % Visualize
                if opts.visualize
                    figure(s); imagesc(seg); axis image off;
                    [yy,xx] = find(pmap > 0);
                    rr = rmap(pmap > 0);
                    assert(numel(xx) == numel(yy) && numel(xx) == numel(rr))
                    viscircles([xx,yy],rr,'Color','r','EnhanceVisibility',true,'Linewidth',0.5);
                end
            end
            progress('Computing MAT...',i,nImages,ticStart,-1);
        end
        BMAX500.(set{1}) = matgt;
    end
    save(savePath, 'BMAX500');
end



