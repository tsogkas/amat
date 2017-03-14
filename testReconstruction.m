function testReconstruction(models,varargin)

% Default testing options ------------------------------------------------
opts = {'dataset',      'BSDS500',...
        'set',          'val',...   % 'val' or 'test'
        'visualize',    false,...
        'parpoolSize',  feature('numcores')  % set to 0 to run serially
       };                        
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths();
if ischar(opts.set) && strcmp(opts.set, 'val')
    opts.imPath    = fullfile(paths.bsds500im,'val');
    opts.gtPath    = fullfile(paths.bsds500gt,'val');
    imageList = dir(fullfile(opts.imPath, '*jpg'));
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    opts.imPath    = fullfile(paths.bsds500im,'test');
    opts.gtPath    = fullfile(paths.bsds500gt,'test');
    imageList = dir(fullfile(opts.imPath, '*jpg'));
elseif isstruct(opts.set)
    disp('Data provided in struct form')
    imageList = opts.set;
    if strcmp(opts.dataset, 'BSDS500')
        if numel(imageList) == 100, opts.set = 'val'; else opts.set = 'test'; end
    end
else
    error('set can be ''val'', ''test'', or a struct containing test data')
end

% Load models and initialize stats ----------------------------------------
if ~iscell(models), models = {models}; end
for m=1:numel(models) 
    models{m} = reconstructDataset(models{m},imageList,opts,paths);
end

% Store stats and save to disk
for m=1:numel(models)
    models{m}.(opts.dataset).(opts.set).stats = models{m}.stats;
    models{m}.(opts.dataset).(opts.set).opts = opts;
    models{m} = rmfield(models{m},'stats');
    modelPath = fullfile(paths.amat.models, models{m}.name);
    model = models{m}; save(modelPath, 'model')
end

% -------------------------------------------------------------------------
function model = reconstructDataset(model,imageList,opts,paths)
% -------------------------------------------------------------------------
switch lower(model)
    case {'amat','gtseg'}
        model = struct('name',model);
    otherwise % load MIL or CNN model
        model = loadModelFromMatFile(model,paths);
end

% Initialize stats
opts.nImages = numel(imageList);
MSE  = zeros(opts.nImages,1);
PSNR = zeros(opts.nImages,1);
SSIM = zeros(opts.nImages,1);
COMP = zeros(opts.nImages,1);

% Evaluate models on test images and compute approximate reconstructions --
modelName = lower(model.name);
ticStart = tic;
parfor (i=1:opts.nImages, opts.parpoolSize)
% for i=1:opts.nImages
    if isfield(imageList(i), 'isdir')
        img = imread(fullfile(opts.imPath,imageList(i).name));
        [~,iid] = fileparts(imageList(i).name);
        gt = load(fullfile(opts.gtPath,[iid '.mat' ])); gt = gt.groundTruth;
        seg  = zeros([size(gt{1}.Segmentation), numel(gt)],'uint16');
        for s=1:numel(gt), seg(:,:,s) = gt{s}.Segmentation; end
    else 
        img = imageList(i).img;
        seg = imageList(i).seg;
    end
    
    % Subsample image for efficiency
    imgResized = imresize(img, 0.5);
    
    switch modelName
        case 'amat'
            [rec,COMP(i)] = reconstructionAMAT(imresize(L0Smoothing(img),0.5));
        case 'gtseg'
            [rec,COMP(i)] = reconstructionGTSEG(imgResized, imresize(seg,0.5,'nearest'));
        case 'deepskel'
            [rec,COMP(i)] = reconstructionDeepSkel(model,imgResized);
        otherwise % MIL
            [rec,COMP(i)] = reconstructionMIL(model,imgResized);
    end
        
    MSE(i)  = immse(double(rec),im2double(imgResized));
    PSNR(i) = psnr(double(rec), im2double(imgResized));
    SSIM(i) = ssim(double(rec), im2double(imgResized));
    msg = sprintf('Testing image reconstruction on %s %s set. ', opts.dataset, opts.set);
    progress(msg,i,opts.nImages,ticStart,-1);
end

% Store stats
model.stats.mse  = MSE;
model.stats.psnr = PSNR;
model.stats.ssim = SSIM;
model.stats.compression = COMP;

% -------------------------------------------------------------------------
function [rec,comp] = reconstructionAMAT(img)
% -------------------------------------------------------------------------
mat = amat(img);
mat.branches = groupMedialPoints(mat);
mat = refineMAT(mat);
rec = mat.reconstruction;
% Sometimes the reconstruction of the simplified AMAT has a few holes, so
% we use inpainting to fill them.
if nnz(isnan(rec))  
    rec = reshape(inpaint_nans(double(rec)), size(rec,1), size(rec,2), []);
    rec = min(1, max(0,rec));
end
comp = nnz(any(mat.axis,3))/(size(img,1)*size(img,2));

% -------------------------------------------------------------------------
function rec = reconstructionDeepSkel(model,img)
% -------------------------------------------------------------------------
rfields = [0,14,40,92,196];
net = model.net;
img = bsxfun(@minus, single(img), reshape(net.meta.averageImage,1,1,[]));
net.eval({'input',img});
lout= net.getLayer(net.getLayerIndex('softmax'));
spb = net.vars(lout.outputIndexes); % symmetry probabilities
spb = spb.value > threshold; % set this to the optimal threshold (contained in model struct)
[~,scales] = max(spb.value,[],3); % estimate scales
scales = rfields(scales);

% -------------------------------------------------------------------------
function [rec,comp] = reconstructionMIL(model,img)
% -------------------------------------------------------------------------
histf = computeHistogramFeatures(img);
spb = spbMIL(img, 'featureSet',model.opts.featureSet,'w',model.w,'histFeatures',histf);
% Create rectangular filters at all scales and orientations
filters = cell(numel(histf.thetas), numel(histf.scales));
for s=1:numel(histf.scales)
    sc = histf.scales(s);
    rect = ones(2*sc+1, 2*histf.opts.ratio*sc+1);
    for o=1:numel(histf.thetas)
        filters{o,s} = imrotate(rect,rad2deg(histf.thetas(o)));
    end
end
% Pad input image, scales and thetas by replicating image border. This is
% necessary to compute mean value encoding correctly, and to avoid out of
% bound errors. Pad = 2*max(scales) to account for rotated versions of filters.
img    = im2double(img);
pad    = round(max(vec(cellfun(@(x)max(size(x)), filters)))/2);
imgPad = padarray(img,[pad,pad],'replicate','both');
pb     = padarray(spb.thin,[pad,pad],0,'both');
scales = padarray(spb.scalesMap,[pad,pad],1,'both');
thetas = padarray(spb.orientMap,[pad,pad],1,'both');
% Sort medial point locations, scales and orientations by their scores.
[pbSorted, indSorted] = sort(pb(:),'descend');
scales = scales(indSorted);
thetas = thetas(indSorted);

% Selecting medial points with the highest scores, compute mean values over
% the respective filter support and stop when the entire image has been covered.
[H,W,C] = size(imgPad);
rec = zeros(H,W,C);
covered = zeros(H,W); covered(border(covered,pad)) = 1;
[y,x] = ind2sub([H,W], indSorted);
i = 1;
while pbSorted(i) > 0.1 && ~all(covered(:))
% while pbSorted(i) > 0.1
    yy = y(i); xx = x(i);
    sc = scales(i);
    th = thetas(i);
    rf = filters{th,sc};
    [h,w] = size(rf);
    hs = floor(h/2); y1 = yy-hs; y2 = y1+h-1;
    ws = floor(w/2); x1 = xx-ws; x2 = x1+w-1;
    assert(y2-y1+1 == h);
    assert(x2-x1+1 == w);
    % Compute mean value encoding of local rectangular patch
    patch = bsxfun(@times, imgPad(y1:y2,x1:x2,:), rf);
    mval = sum(sum(patch))/nnz(rf);
    % Only fill in pixels that have no<t been already covered
    patch = bsxfun(@times, mval, rf );
    rec(y1:y2,x1:x2,:) = rec(y1:y2,x1:x2,:)+patch;
    covered(y1:y2,x1:x2) = covered(y1:y2,x1:x2) + rf;
    i = i+1;
end
rec = bsxfun(@rdivide, rec, covered);
rec = rec(pad+1:end-pad,pad+1:end-pad,:);
% Fill in holes that might be left
rec = reshape(inpaint_nans(rec), size(img,1), size(img,2),3);
comp = i / (size(img,1)*size(img,2));

% -------------------------------------------------------------------------
function [rec,comp] = reconstructionGTSEG(img,seg)
% -------------------------------------------------------------------------
% Baseline reconstruction statistics using "oracle" segmentations from BSDS
img       = im2double(img);
numSegs   = size(seg,3);
nnzPixels = zeros(numSegs,1);
rec       = zeros([size(img), numSegs]);
for i=1:numSegs
    segi = seg(:,:,i);
    labels = unique(segi);
    numSegments = numel(labels);
    for s=1:numSegments
        segment = segi == labels(s);
        area = nnz(segment);
        boundaries = bwperim(segment);
        nnzPixels(i) = nnzPixels(i) + nnz(boundaries);
        meanColor = sum(sum(bsxfun(@times, img, double(segment))))/area;
        meanSegment = bsxfun(@times, double(repmat(segment,1,1,3)), meanColor);
        rec(:,:,:,i) = rec(:,:,:,i) + meanSegment; 
    end
end
% Select segmentation with maximum compression and return respective rec
[comp,idx] = min(nnzPixels / (size(img,1)*size(img,2)));
rec = rec(:,:,:,idx);



% -------------------------------------------------------------------------
function model = loadModelFromMatFile(model,paths)
% -------------------------------------------------------------------------
if exist(fullfile(paths.spbmil.models, model),'file') % MIL model
    tmp = load(fullfile(paths.spbmil.models, model));
elseif exist(model,'file')
    tmp = load(model);
end
if isfield(tmp,'model')     % MIL model
    model = tmp.model;
elseif isfield(tmp,'net')   % DagNN model
    model = struct('trainStats',tmp.stats, 'name','deepskel',...
                   'net',dagnn.DagNN.loadobj(tmp.net));
end
