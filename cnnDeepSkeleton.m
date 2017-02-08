function [net, info] = cnnDeepSkeleton(varargin)

% Network options
opts.dataDir = [];
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
opts.expDir = fullfile(paths.deepskel.output, ['deepskel-' sfx]) ;
[opts, varargin] = vl_argparse(opts, varargin) ;

% GPU options
opts.numFetchThreads = 12 ;
opts.train = struct() ;
opts = vl_argparse(opts, varargin) ;
if ~isfield(opts.train, 'gpus'), opts.train.gpus = []; end;

% Prepare training data
imdb = getBMAX500Imdb('dataDir', opts.dataDir);

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
modelPath = fullfile(opts.expDir, 'net-deployed.mat')

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

% Images in BSDS500 (hence in BMAX500 too) have dimensions that are either 
% 321x481 or 481x321. Accounting for the scales and orientations used in
% data augmentation, the maximum dimensions for the input data are
% ceil(481*1.2) = 578. For simplicity and speed while getting batches, we
% use a fixed size for the input images and label maps, by centering the
% (possibly smaller) image in this 578x578 space.
maxH = ceil(max(scales)*481); maxW = maxH;

% Load images/groundtruth and augment scales
% For each set
sets = {'train','val','test'};
for s=1:numel(sets)
    set = bmax500.(sets{s});
    % First count the total number of groundtruth maps for preallocation
    nImages = 0;
    for i=1:numel(set)
        nImages = nImages + size(set(i).pts,3);
    end
    % Preallocate
    imgs = zeros(maxH,maxW,3,nImages*(1+numel(scales)),'uint8');
    radi = zeros(maxH,maxW,1,nImages*(1+numel(scales)),'uint8');
    lbls = false(maxH,maxW,1,nImages*(1+numel(scales)));
    % Now assemble all examples
    ind = 1;
    for i=1:numel(set)
        ngt = size(set(i).pts,3); % #gt for each image
        for scale = scales  % create scaled versions of images and gt
            img = set(i).img;
            pts = set(i).pts;
            rad = set(i).rad;
            if scale ~= 1
                img = imresize(img,scale,'bilinear');
                rad = imresize(rad,scale,'bilinear')*scale;
                pts = imresizeCrisp(pts,scale);
            end
            % Center data in the maxH x maxW array
            [H,W,~] = size(img);
            hs = ceil((maxH-H)/2); ws = ceil((maxW-W)/2);
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

% Augment orientations
for s=1:numel(sets)
    imgs = imdb.(sets{s}).images;
    lbls = imdb.(sets{s}).labels;
    radi = imdb.(sets{s}).radius;
    for theta = thetas
        if theta ~= 0
            imdb.(sets{s}).images = cat(4,imdb.(sets{s}).images,imrotate(imgs,theta));
            imdb.(sets{s}).labels = cat(4,imdb.(sets{s}).labels,imrotate(lbls,theta));
            imdb.(sets{s}).radius = cat(4,imdb.(sets{s}).radius,imrotate(radi,theta));
        end
    end
end
clear imgs lbls radi

% Tidy imdb
nTrain = size(imdb.train.images,4);
nVal   = size(imdb.val.images,4);
nTest  = size(imdb.tesst.images,4);
imdb.images.data = cat(4, imdb.train.images, imdb.val.images, imdb.test.images);
imdb.train.images = []; imdb.val.images = []; imdb.test.images = [];
imdb.images.labels = cat(4, imdb.train.labels, imdb.val.labels, imdb.test.labels);
imdb.train.labels = []; imdb.val.labels = []; imdb.test.labels = [];
imdb.images.radius = cat(4, imdb.train.radius, imdb.val.radius, imdb.test.radius);
imdb.train.radius = []; imdb.val.radius = []; imdb.test.radius = [];
imdb = rmfield(imdb, {'train','val','test'});
imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal) 3*ones(1,nTest)];
imdb.meta.sets = sets;

% Subtract mean image
imdb.images.data = bsxfun(@minus, imdb.images.data, []);









