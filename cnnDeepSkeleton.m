function [net, info] = cnnDeepSkeleton(varargin)

% Network options
opts.mode = 'train';
opts.dataDir = constructBMAX500('resize',0.5);
opts.modelType = 'vgg-vd-16' ;
opts.network = [] ;
opts.networkType = 'dagnn' ;
opts.batchNormalization = true ;
opts.weightInitMethod = 'gaussian' ;
[opts, varargin] = vl_argparse(opts, varargin) ;

% Suffix and results directory setup
paths = setPaths();
sfx = opts.modelType ;
if opts.batchNormalization, sfx = [sfx '-bnorm'] ; end
sfx = [sfx '-' opts.networkType] ;
opts.expDir = fullfile(paths.amat.output, ['deepskel-' sfx]) ;
[opts, varargin] = vl_argparse(opts, varargin) ;

% GPU options
opts.numFetchThreads = 12 ;
opts.train = struct() ;
opts = vl_argparse(opts, varargin) ;
if ~isfield(opts.train, 'gpus'), opts.train.gpus = []; end;

% Prepare training data
imdb = getBMAX500Imdb('dataDir', opts.dataDir, 'mode',opts.mode);

% Initialize model
net = cnnInit('model', opts.modelType, ...
              'batchNormalization', opts.batchNormalization, ...
              'weightInitMethod', opts.weightInitMethod, ...
              'networkType', opts.networkType, ...
              'averageImage', rgbMean, ...
              'colorDeviation', rgbDeviation, ...
              'classNames', imdb.classes.name, ...
              'classDescriptions', imdb.classes.description) ;


% Train          
trainFn = @cnn_train_dag ;
[net, info] = trainFn(net, imdb, getBatchFn(opts, net.meta), ...
                      'expDir', opts.expDir, ...
                      net.meta.trainOpts, ...
                      opts.train) ;

% -------------------------------------------------------------------------
%                                                                    Deploy
% -------------------------------------------------------------------------

net = cnn_imagenet_deploy(net) ;
modelPath = fullfile(opts.expDir, 'net-deployed.mat');

switch opts.networkType
  case 'simplenn'
    save(modelPath, '-struct', 'net') ;
  case 'dagnn'
    net_ = net.saveobj() ;
    save(modelPath, '-struct', 'net_') ;
    clear net_ ;
end

% -------------------------------------------------------------------------
function fn = getBatchFn(opts, meta)
% -------------------------------------------------------------------------

if numel(meta.normalization.averageImage) == 3
  mu = double(meta.normalization.averageImage(:)) ;
else
  mu = imresize(single(meta.normalization.averageImage), ...
                meta.normalization.imageSize(1:2)) ;
end

useGpu = numel(opts.train.gpus) > 0 ;

bopts.test = struct(...
  'useGpu', useGpu, ...
  'numThreads', opts.numFetchThreads, ...
  'imageSize',  meta.normalization.imageSize(1:2), ...
  'cropSize', meta.normalization.cropSize, ...
  'subtractAverage', mu) ;

% Copy the parameters for data augmentation
bopts.train = bopts.test ;
for f = fieldnames(meta.augmentation)'
  f = char(f) ;
  bopts.train.(f) = meta.augmentation.(f) ;
end

fn = @(x,y) getBatch(bopts,useGpu,lower(opts.networkType),x,y) ;

% -------------------------------------------------------------------------
function varargout = getBatch(opts, useGpu, networkType, imdb, batch)
% -------------------------------------------------------------------------
images = strcat([imdb.imageDir filesep], imdb.images.name(batch)) ;
if ~isempty(batch) && imdb.images.set(batch(1)) == 1
  phase = 'train' ;
else
  phase = 'test' ;
end
data = getImageBatch(images, opts.(phase), 'prefetch', nargout == 0) ;
if nargout > 0
  labels = imdb.images.label(batch) ;
  switch networkType
    case 'simplenn'
      varargout = {data, labels} ;
    case 'dagnn'
      varargout{1} = {'input', data, 'label', labels} ;
  end
end

% -------------------------------------------------------------------------
function imdb = getBMAX500Imdb(varargin)
% -------------------------------------------------------------------------
opts.mode = 'train';
opts.dataDir = [] ;
opts.lite = false ;
opts = vl_argparse(opts, varargin) ;
if ischar(opts.dataDir) && exist(opts.dataDir,'file')
    tmp = load(opts.dataDir); bmax500 = tmp.BMAX500; clear tmp;
elseif isstruct(opts.dataDir)
    bmax500 = opts.dataDir;
else
    error('Input should be a path to the data or a struct containing the data')
end

% Data augmentation parameters
scales = [0.8 1 1.2];    % remember to scale the radii accordingly!
thetas = [0,90,180,270];  % angles in degrees

% For simplicity and speed while getting batches, we use a fixed size for 
% the input images and label maps, by centering the (possibly smaller) 
% image in a D x D space. D is computed based on the maximum image
% dimensions in the dataset and the maximum scaling used.
% NOTE: Checking one image is enough because all images are either DxK or
% KxD.
Dmax = max(size(bmax500.train(1).img)); 
Hmax = ceil(max(scales)*Dmax); Wmax = Hmax;

% Load only the necessary subsets to save RAM
if strcmp(opts.mode, 'test')
    sets = {'test'};
else
    sets = {'train','val'}; 
end

% Load images/groundtruth and augment scales
for s=1:numel(sets)
    set = bmax500.(sets{s});
    % First count the total number of groundtruth maps for preallocation
    nImages = 0;
    for i=1:numel(set)
        nImages = nImages + size(set(i).pts,3);
    end
    % Preallocate
    imgs = zeros(Hmax,Wmax,3,nImages*(1+numel(scales)),'uint8');
    radi = zeros(Hmax,Wmax,1,nImages*(1+numel(scales)),'uint8');
    lbls = false(Hmax,Wmax,1,nImages*(1+numel(scales)));
    % Now assemble all examples
    ind = 1;
    for i=1:numel(set)
        ngt = size(set(i).pts,3); % #gt for each image
        for scale = scales  % create scaled versions of images and gt
            img = set(i).img;
            pts = set(i).pts;
            rad = set(i).rad;
            if scale ~= 1 && (strcmp(sets{s},'train') || strcmp(sets{s},'trainval'))
                img = imresize(img,scale,'bilinear');
                rad = imresize(rad,scale,'bilinear')*scale;
                pts = imresizeCrisp(pts,scale);
            end
            % Center data in the maxH x maxW array
            [H,W,~] = size(img);
            hs = ceil((Hmax-H+1)/2); ws = ceil((Wmax-W+1)/2);
            imgs(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = repmat(img,[1,1,1,ngt]);
            radi(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = reshape(rad,H,W,1,[]);
            lbls(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = reshape(pts,H,W,1,[]);
            ind = ind+ngt;
        end
    end
    imdb.(sets{s}).images = imgs;
    imdb.(sets{s}).labels = lbls;
    imdb.(sets{s}).radius = radi;
end

% Tidy imdb
switch opts.mode
    case 'test'
        imdb.images.data   = imdb.test.images;  imdb.test.images = [];
        imdb.images.labels = imdb.test.labels;  imdb.test.labels = [];
        imdb.images.radius = imdb.test.radius;  imdb.test.radius = [];
        imdb = rmfield(imdb, 'test');
    case 'train'
        imdb.images.data   = imdb.train.images; imdb.train.images = [];
        imdb.images.labels = imdb.train.labels; imdb.train.labels = []; 
        imdb.images.radius = imdb.train.radius; imdb.train.radius = []; 
        imdb = rmfield(imdb, 'train');
    case 'trainval'
        imdb.images.data = cat(4, imdb.train.images, imdb.val.images);
        imdb.train.images = []; imdb.val.images = [];
        imdb.images.labels = cat(4, imdb.train.labels, imdb.val.labels);
        imdb.train.labels = []; imdb.val.labels = [];
        imdb.images.radius = cat(4, imdb.train.radius, imdb.val.radius);
        imdb.train.radius = []; imdb.val.radius = [];
        imdb = rmfield(imdb, {'train','val'});
    otherwise, error('opts.mode can be ''train'',''trainval'', or''test''')
end

% Augment orientations
imgs = imdb.images.data;
lbls = imdb.images.labels;
radi = imdb.images.radius;
for theta = thetas
    if theta ~= 0
        imdb.images.data   = cat(4, imdb.images.data,   imrotate(imgs,theta));
        imdb.images.labels = cat(4, imdb.images.labels, imrotate(lbls,theta));
        imdb.images.radius = cat(4, imdb.images.radius, imrotate(radi,theta));
    end
end
clear imgs lbls radi img pts rad

% Augment flipping
imdb.images.data = cat(4, imdb.images.data,...
    fliplr(imdb.images.data), flipud(imdb.images.data));
imdb.images.labels = cat(4, imdb.images.labels,...
    fliplr(imdb.images.labels), flipud(imdb.images.labels));
imdb.images.radius = cat(4, imdb.images.radius,...
    fliplr(imdb.images.radius), flipud(imdb.images.radius));

% Assign set indexes
if strcmp(opts.mode, 'train')
    nTrain = size(imdb.images.data,4);
    nVal = size(imdb.val.images,4);
    imdb.images.data = cat(4, imdb.images.data, imdb.val.images);
    imdb.val.images = []; imdb = rmfield(imdb,'val');
    imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal)];
elseif strcmp(opts.mode, 'trainval') % what to do for validation here??
    nTrain = size(imdb.images.data,4);
    imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal)];
else error('Invalid opts.mode')
end
imdb.meta.sets = sets;

% Subtract mean image
imdb.images.data = bsxfun(@minus, imdb.images.data, []);

% -------------------------------------------------------------------------
function net = cnnInit(varargin)
% -------------------------------------------------------------------------
opts.scale = 1 ;
opts.initBias = 0 ;
opts.weightDecay = 1 ;
opts.weightInitMethod = 'gaussian' ;
opts.model = 'alexnet' ;
opts.batchNormalization = false ;
opts.networkType = 'simplenn' ;
opts.cudnnWorkspaceLimit = 1024*1024*1204 ; % 1GB
opts.classNames = {} ;
opts.classDescriptions = {} ;
opts.averageImage = zeros(3,1) ;
opts.colorDeviation = zeros(3) ;
opts = vl_argparse(opts, varargin) ;

net.meta.normalization.imageSize = [224, 224, 3] ;
net = vgg_vd(net, opts) ;
bs = 32 ;

% final touches
switch lower(opts.weightInitMethod)
  case {'xavier', 'xavierimproved'}
    net.layers{end}.weights{1} = net.layers{end}.weights{1} / 10 ;
end
net.layers{end+1} = struct('type', 'softmaxloss', 'name', 'loss') ;

% Meta parameters
net.meta.inputSize = [net.meta.normalization.imageSize, 32] ;
net.meta.normalization.cropSize = net.meta.normalization.imageSize(1) / 256 ;
net.meta.normalization.averageImage = opts.averageImage ;
net.meta.classes.name = opts.classNames ;
net.meta.classes.description = opts.classDescriptions;
net.meta.augmentation.jitterLocation = true ;
net.meta.augmentation.jitterFlip = true ;
net.meta.augmentation.jitterBrightness = double(0.1 * opts.colorDeviation) ;
net.meta.augmentation.jitterAspect = [2/3, 3/2] ;

if ~opts.batchNormalization
  lr = logspace(-2, -4, 60) ;
else
  lr = logspace(-1, -4, 20) ;
end

net.meta.trainOpts.learningRate = lr ;
net.meta.trainOpts.numEpochs = numel(lr) ;
net.meta.trainOpts.batchSize = bs ;
net.meta.trainOpts.weightDecay = 0.0005 ;

% Fill in default values
net = vl_simplenn_tidy(net) ;

% Switch to DagNN if requested
switch lower(opts.networkType)
  case 'simplenn'
    % done
  case 'dagnn'
    net = dagnn.DagNN.fromSimpleNN(net, 'canonicalNames', true) ;
    net.addLayer('top1err', dagnn.Loss('loss', 'classerror'), ...
                 {'prediction','label'}, 'top1err') ;
    net.addLayer('top5err', dagnn.Loss('loss', 'topkerror', ...
                                       'opts', {'topK',5}), ...
                 {'prediction','label'}, 'top5err') ;
  otherwise
    assert(false) ;
end

% --------------------------------------------------------------------
function net = add_block(net, opts, id, h, w, in, out, stride, pad)
% --------------------------------------------------------------------
info = vl_simplenn_display(net) ;
fc = (h == info.dataSize(1,end) && w == info.dataSize(2,end)) ;
if fc
  name = 'fc' ;
else
  name = 'conv' ;
end
convOpts = {'CudnnWorkspaceLimit', opts.cudnnWorkspaceLimit} ;
net.layers{end+1} = struct('type', 'conv', 'name', sprintf('%s%s', name, id), ...
                           'weights', {{init_weight(opts, h, w, in, out, 'single'), ...
                             ones(out, 1, 'single')*opts.initBias}}, ...
                           'stride', stride, ...
                           'pad', pad, ...
                           'dilate', 1, ...
                           'learningRate', [1 2], ...
                           'weightDecay', [opts.weightDecay 0], ...
                           'opts', {convOpts}) ;
if opts.batchNormalization
  net.layers{end+1} = struct('type', 'bnorm', 'name', sprintf('bn%s',id), ...
                             'weights', {{ones(out, 1, 'single'), zeros(out, 1, 'single'), ...
                               zeros(out, 2, 'single')}}, ...
                             'epsilon', 1e-4, ...
                             'learningRate', [2 1 0.1], ...
                             'weightDecay', [0 0]) ;
end
net.layers{end+1} = struct('type', 'relu', 'name', sprintf('relu%s',id)) ;

% -------------------------------------------------------------------------
function weights = init_weight(opts, h, w, in, out, type)
% -------------------------------------------------------------------------
% See K. He, X. Zhang, S. Ren, and J. Sun. Delving deep into
% rectifiers: Surpassing human-level performance on imagenet
% classification. CoRR, (arXiv:1502.01852v1), 2015.

switch lower(opts.weightInitMethod)
  case 'gaussian'
    sc = 0.01/opts.scale ;
    weights = randn(h, w, in, out, type)*sc;
  case 'xavier'
    sc = sqrt(3/(h*w*in)) ;
    weights = (rand(h, w, in, out, type)*2 - 1)*sc ;
  case 'xavierimproved'
    sc = sqrt(2/(h*w*out)) ;
    weights = randn(h, w, in, out, type)*sc ;
  otherwise
    error('Unknown weight initialization method''%s''', opts.weightInitMethod) ;
end

% --------------------------------------------------------------------
function net = add_dropout(net, opts, id)
% --------------------------------------------------------------------
if ~opts.batchNormalization
  net.layers{end+1} = struct('type', 'dropout', ...
                             'name', sprintf('dropout%s', id), ...
                             'rate', 0.5) ;
end

% --------------------------------------------------------------------
function net = vgg_vd(net, opts)
% --------------------------------------------------------------------

net = 
net.layers = {} ;
net = add_block(net, opts, '1_1', 3, 3, 3, 64, 1, 1) ;
net = add_block(net, opts, '1_2', 3, 3, 64, 64, 1, 1) ;
net.layers{end+1} = struct('type', 'pool', 'name', 'pool1', ...
                           'method', 'max', ...
                           'pool', [2 2], ...
                           'stride', 2, ...
                           'pad', 0) ;

net = add_block(net, opts, '2_1', 3, 3, 64, 128, 1, 1) ;
net = add_block(net, opts, '2_2', 3, 3, 128, 128, 1, 1) ;
net.layers{end+1} = struct('type', 'pool', 'name', 'pool2', ...
                           'method', 'max', ...
                           'pool', [2 2], ...
                           'stride', 2, ...
                           'pad', 0) ;

net = add_block(net, opts, '3_1', 3, 3, 128, 256, 1, 1) ;
net = add_block(net, opts, '3_2', 3, 3, 256, 256, 1, 1) ;
net = add_block(net, opts, '3_3', 3, 3, 256, 256, 1, 1) ;
net.layers{end+1} = struct('type', 'pool', 'name', 'pool3', ...
                           'method', 'max', ...
                           'pool', [2 2], ...
                           'stride', 2, ...
                           'pad', 0) ;

net = add_block(net, opts, '4_1', 3, 3, 256, 512, 1, 1) ;
net = add_block(net, opts, '4_2', 3, 3, 512, 512, 1, 1) ;
net = add_block(net, opts, '4_3', 3, 3, 512, 512, 1, 1) ;
net.layers{end+1} = struct('type', 'pool', 'name', 'pool4', ...
                           'method', 'max', ...
                           'pool', [2 2], ...
                           'stride', 2, ...
                           'pad', 0) ;

net = add_block(net, opts, '5_1', 3, 3, 512, 512, 1, 1) ;
net = add_block(net, opts, '5_2', 3, 3, 512, 512, 1, 1) ;
net = add_block(net, opts, '5_3', 3, 3, 512, 512, 1, 1) ;
net.layers{end+1} = struct('type', 'pool', 'name', 'pool5', ...
                           'method', 'max', ...
                           'pool', [2 2], ...
                           'stride', 2, ...
                           'pad', 0) ;

net = add_block(net, opts, '6', 7, 7, 512, 4096, 1, 0) ;
net = add_dropout(net, opts, '6') ;

net = add_block(net, opts, '7', 1, 1, 4096, 4096, 1, 0) ;
net = add_dropout(net, opts, '7') ;

net = add_block(net, opts, '8', 1, 1, 4096, 1000, 1, 0) ;
net.layers(end) = [] ;
if opts.batchNormalization, net.layers(end) = [] ; end