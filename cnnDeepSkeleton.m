function [net, info] = cnnDeepSkeleton(varargin)

% Network options
opts.mode = 'train';
opts.dataDir = constructBMAX500('resize',0.5);
opts.modelType = 'vgg-vd-16' ;
opts.network = [] ;
opts.networkType = 'dagnn' ;
opts.batchNormalization = true ;
opts.weightInitMethod = 'gaussian' ;
opts.averageImage = [116.66877; 122.67892; 104.00699];  
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

% Initialize model
net = cnnInit('initNetPath', paths.vgg16);

% Prepare training data
imdb = getBMAX500Imdb('dataDir', opts.dataDir,...
                      'mode',opts.mode,...
                      'averageImage', opts.averageImage);

% Train          
[net, info] = cnn_train_dag(net, imdb, @getBatch, ...
                      'expDir', opts.expDir, ...
                      net.meta.trainOpts, ...
                      opts.train) ;


% -------------------------------------------------------------------------
function out = getBatch(imdb, batch)
% -------------------------------------------------------------------------
images = single(imdb.images.data(batch));
labels = single(imdb.images.label(batch));
images = bsxfun(@minus, images, reshape(imdb.meta.averageImage,1,1,3));
out{1} = {'input', images, 'label', labels} ;

% -------------------------------------------------------------------------
function imdb = getBMAX500Imdb(varargin)
% -------------------------------------------------------------------------
opts.mode = 'train';
opts.dataDir = [] ;
opts.lite = false ;
opts.averageImage = zeros(3,1);
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

% Load images/groundtruth and augment scales. Each segmentation is
% considered as an individual training example. 
% TODO: is there a better way to handle multiple segmentations?
for s=1:numel(sets)
    set = bmax500.(sets{s});
    % First count the total number of groundtruth maps for preallocation
    nImages = sum(cellfun(@(x) size(x,3), {set(:).pts}));
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
            for k=1:size(pts,3)
                plotDisks(img,pts(:,:,k),rad(:,:,k));
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
        imdb.images.data   = cat(4, imdb.images.data,   imrotate(imgs,theta,'bilinear'));
        imdb.images.labels = cat(4, imdb.images.labels, imrotate(lbls,theta,'nearest'));
        imdb.images.radius = cat(4, imdb.images.radius, imrotate(radi,theta,'nearest'));
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
elseif strcmp(opts.mode, 'trainval') % TODO: what to do for validation here??
    nTrain = size(imdb.images.data,4);
    imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal)];
else error('Invalid opts.mode')
end
imdb.meta.sets = sets;
imdb.meta.averageImage = opts.averageImage;

% -------------------------------------------------------------------------
function net = cnnInit(varargin)
% -------------------------------------------------------------------------
% TODO: We want to fine-tune, starting from the vgg-16 initialization so we
% have to:
% 1) Load the weights of that network (download it from the matconvet website)
% 2) Set the values for the learning parameters to the same values as the
% ones used in the DeepSkeleton paper
% (https://github.com/zeakey/DeepSkeleton/blob/master/examples/DeepSkeleton/train_val.prototxt)

opts.initNetPath = [];
opts.cudnnWorkspaceLimit = 1024*1024*1204 ; % 1GB
opts = vl_argparse(opts, varargin) ;

% Define network architecture and initialize parameters with random values
net = vgg16deepskel();
net.initParams(); 

% Set learning parameters for Convolutional layers
vgg = load(opts.initNetPath); 
for l=1:numel(vgg.layers)
    layer = vgg.layers{l};
    if strcmp(layer.type, 'conv')
        % Initialize weights from the vgg-16 network
        indw = net.getParamIndex([layer.name '_w']);
        indb = net.getParamIndex([layer.name '_b']);
        if isnan(indw) || isnan(indb)
            warning(['Skipping non-existent layer ' layer.name '...'])
        else
            % Set learning parameters according to the deep skeleton net.
            net.params(indw).value = layer.weights{1};
            net.params(indb).value = layer.weights{2};
            net.params(indw).weightDecay  = 1;
            net.params(indb).weightDecay  = 0;
            if ismember(layer.name, {'conv5_1','conv5_2','conv5_3'})
                net.params(indw).learningRate = 100;
                net.params(indb).learningRate = 200;
            else
                net.params(indw).learningRate = 1;
                net.params(indb).learningRate = 2;
            end
        end
    end
end

% Learning parameters for deconvolutional (upsampling) layers
for l={'deconv2','deconv3','deconv4','deconv5'}
    indw = net.getParamIndex([l{1} '_w']);
    indb = net.getParamIndex([l{1} '_b']);
    net.params(indw).learningRate = 0;
    net.params(indw).weightDecay  = 1;
    net.params(indb).learningRate = 0;
    net.params(indb).weightDecay  = 0;
end

% Learning parameters for side output layers
for l={'score_dsn2','score_dsn3','score_dsn4','score_dsn5'}
    indw = net.getParamIndex([l{1} '_w']);
    indb = net.getParamIndex([l{1} '_b']);
    net.params(indw).learningRate = 0.01;
    net.params(indw).weightDecay  = 1;
    net.params(indb).learningRate = 0.02;
    net.params(indb).weightDecay  = 0;
end

% Layers for score fusion layers
for l={'cat0_score','cat1_score','cat2_score','cat3_score','cat4_score'}
    indw = net.getParamIndex([l{1} '_w']);
    indb = net.getParamIndex([l{1} '_b']);
    net.params(indw).learningRate = 0.05; % 0.01 for cat2-score in the original code - probably a typo
    net.params(indw).weightDecay  = 1;
    net.params(indb).learningRate = 0.002;
    net.params(indb).weightDecay  = 0;
    % Initialize score fusion to uniform weighting
    switch l{1}
        case {'cat0_score','cat1_score'}
            val = 0.25;
        case 'cat2_score'
            val = 1/3;
        case 'cat3_score'
            val = 0.5;
        case 'cat4_score'
            val = 1;
    end
    net.params(indw).value = val * ones(size(net.params(indw).value));
end

% Meta parameters 
% - BSDS trainval set: 300 (200 for train set)
% - Augmentation factor: 36 (~36*6 = 216 if we consider each segmentation 
%   as a separate training example).
% - Batchsize: 10
% Hence: #iters/epoch = 300*36 / 10 = 1080, and we need ~15-20 epochs to
% reach 20K iterations (~3 epochs in the case we use all segmentations).
net.meta.trainOpts.learningRate = [1e-6*ones(1,5), 1e-7*ones(1,5), 1e-8*ones(1,5)] ;
net.meta.trainOpts.numEpochs = numel(net.meta.trainOpts.learningRate) ;
net.meta.trainOpts.momentum = 0.9 ;
net.meta.trainOpts.weightDecay = 0.0002 ;
net.meta.trainOpts.batchSize = 10 ;

% --------------------------------------------------------------------
function net = vgg16deepskel()
% --------------------------------------------------------------------
% Original HED architecture (https://arxiv.org/abs/1504.06375), extended
% with deep side, loss, and fusion layers for combining skeleton responses
% at different scales.

net = dagnn.DagNN();

% Conv Layers -------------------------------------------------------------
% Conv 1 
net.addLayer('conv1_1', dagnn.Conv('size',[3,3,3,64],'stride',1,'pad',1),...
            {'input'}, {'conv1_1'}, {'conv1_1_w', 'conv1_1_b'});
net.addLayer('relu1_1', dagnn.ReLU(), {'conv1_1'}, {'relu1_1'});
net.addLayer('conv1_2', dagnn.Conv('size',[3,3,64,64],'stride',1,'pad',1),...
            {'relu1_1'}, {'conv1_2'}, {'conv1_2_w', 'conv1_2_b'});
net.addLayer('relu1_2', dagnn.ReLU(), {'conv1_2'}, {'relu1_2'});
net.addLayer('pool1', dagnn.Pooling('poolSize',[2 2],'stride',2,'pad',0,...
             'method','max'), {'relu1_2'}, {'pool1'});
% Conv 2 
net.addLayer('conv2_1', dagnn.Conv('size',[3,3,64,128],'stride',1,'pad',1),...
            {'pool1'}, {'conv2_1'},{'conv2_1_w', 'conv2_1_b'});
net.addLayer('relu2_1', dagnn.ReLU(), {'conv2_1'}, {'relu2_1'});        
net.addLayer('conv2_2', dagnn.Conv('size',[3,3,128,128],'stride',1,'pad',1),...
            {'relu2_1'}, {'conv2_2'},{'conv2_2_w', 'conv2_2_b'});
net.addLayer('relu2_2', dagnn.ReLU(), {'conv2_2'}, {'relu2_2'});                
net.addLayer('pool2', dagnn.Pooling('poolSize',[2 2],'stride',2,'pad',0,...
             'method','max') ,{'relu2_2'}, {'pool2'});
                  
% Conv 3 
net.addLayer('conv3_1', dagnn.Conv('size',[3,3,128,256],'stride',1,'pad',1),...
            {'pool2'}, {'conv3_1'},{'conv3_1_w', 'conv3_1_b'});
net.addLayer('relu3_1', dagnn.ReLU(), {'conv3_1'}, {'relu3_1'});
net.addLayer('conv3_2', dagnn.Conv('size',[3,3,256,256],'stride',1,'pad',1),...
            {'relu3_1'}, {'conv3_2'},{'conv3_2_w', 'conv3_2_b'});
net.addLayer('relu3_2', dagnn.ReLU(), {'conv3_2'}, {'relu3_2'});                                
net.addLayer('conv3_3', dagnn.Conv('size',[3,3,256,256],'stride',1,'pad',1),...
            {'relu3_2'}, {'conv3_3'},{'conv3_3_w', 'conv3_3_b'});
net.addLayer('relu3_3', dagnn.ReLU(), {'conv3_2'}, {'relu3_3'});
net.addLayer('pool3', dagnn.Pooling('poolSize',[2 2],'stride',2,'pad',0,...
             'method','max') ,{'relu3_3'}, {'pool3'});
         
% Conv 4 
net.addLayer('conv4_1', dagnn.Conv('size',[3,3,256,512],'stride',1,'pad',1),...
            {'pool3'}, {'conv4_1'},{'conv4_1_w', 'conv4_1_b'});
net.addLayer('relu4_1', dagnn.ReLU(), {'conv4_1'}, {'relu4_1'});
net.addLayer('conv4_2', dagnn.Conv('size',[3,3,512,512],'stride',1,'pad',1),...
            {'relu4_1'}, {'conv4_2'},{'conv4_2_w', 'conv4_2_b'});
net.addLayer('relu4_2', dagnn.ReLU(), {'conv4_2'}, {'relu4_2'});        
net.addLayer('conv4_3', dagnn.Conv('size',[3,3,512,512],'stride',1,'pad',1),...
            {'relu4_2'}, {'conv4_3'},{'conv4_3_w', 'conv4_3_b'});
net.addLayer('relu4_3', dagnn.ReLU(), {'conv4_3'}, {'relu4_3'});                
net.addLayer('pool4', dagnn.Pooling('poolSize',[2 2],'stride',2,'pad',0,...
             'method','max') ,{'relu4_3'}, {'pool4'});
                  
% Conv 5
net.addLayer('conv5_1', dagnn.Conv('size',[3,3,512,512],'stride',1,'pad',1),...
            {'pool4'}, {'conv5_1'},{'conv5_1_w', 'conv5_1_b'});
net.addLayer('relu5_1', dagnn.ReLU(), {'conv5_1'}, {'relu5_1'});                
net.addLayer('conv5_2', dagnn.Conv('size',[3,3,512,512],'stride',1,'pad',1),...
            {'relu5_1'}, {'conv5_2'},{'conv5_2_w', 'conv5_2_b'});
net.addLayer('relu5_2', dagnn.ReLU(), {'conv5_2'}, {'relu5_2'});                
net.addLayer('conv5_3', dagnn.Conv('size',[3,3,512,512],'stride',1,'pad',1),...
            {'relu5_2'}, {'conv5_3'},{'conv5_3_w', 'conv5_3_b'});
net.addLayer('relu5_3', dagnn.ReLU(), {'conv5_3'}, {'relu5_3'});                
net.addLayer('pool5', dagnn.Pooling('poolSize',[2 2],'stride',2,'pad',0,...
             'method','max') ,{'relu5_3'}, {'pool5'});
         
% DSN Layers --------------------------------------------------------------
% DSN 2
net.addLayer('score_dsn2', dagnn.Conv('size',[1,1,128,2]),...
            {'conv2_2'}, {'score_dsn2'},{'score_dsn2_w', 'score_dsn2_b'});
net.addLayer('deconv2', dagnn.ConvTranspose('size',[4,4,2,2],'upsample',2),...
            {'score_dsn2'}, {'score_dsn2_up'}, {'deconv2_w','deconv2_b'});
net.addLayer('loss2', dagnn.Loss(),{'score_dsn2_up','labels'}, {'loss_dsn2'});
                  
% DSN 3
net.addLayer('score_dsn3', dagnn.Conv('size',[1,1,256,3]),...
            {'conv3_3'}, {'score_dsn3'},{'score_dsn3_w', 'score_dsn3_b'});
net.addLayer('deconv3', dagnn.ConvTranspose('size',[8,8,3,3],'upsample',4),...
            {'score_dsn3'}, {'score_dsn3_up'}, {'deconv3_w','deconv3_b'});
net.addLayer('loss3', dagnn.Loss(),{'score_dsn3_up','labels'}, {'loss_dsn3'});
                  
% DSN 4
net.addLayer('score_dsn4', dagnn.Conv('size',[1,1,512,4]),...
            {'conv4_3'}, {'score_dsn4'},{'score_dsn4_w', 'score_dsn4_b'});
net.addLayer('deconv4', dagnn.ConvTranspose('size',[16,16,4,4],'upsample',8),...
            {'score_dsn4'}, {'score_dsn4_up'}, {'deconv4_w','deconv4_b'});
net.addLayer('loss4', dagnn.Loss(),{'score_dsn4_up','labels'}, {'loss_dsn4'});
         
% DSN 5
net.addLayer('score_dsn5', dagnn.Conv('size',[1,1,512,5]),...
            {'conv5_3'}, {'score_dsn5'},{'score_dsn5_w', 'score_dsn5_b'});
net.addLayer('deconv5', dagnn.ConvTranspose('size',[32,32,5,5],'upsample',16),...
            {'score_dsn5'}, {'score_dsn5_up'}, {'deconv5_w','deconv5_b'});
net.addLayer('loss5', dagnn.Loss(),{'score_dsn5_up','labels'}, {'loss_dsn5'});

% Slice side outputs ------------------------------------------------------
net.addLayer('slice2', dagnn.Slice(), {'upscore_dsn2'},...
            {'slice2_0','slice2_1'});
net.addLayer('slice3', dagnn.Slice(), {'upscore_dsn3'},...
            {'slice3_0','slice3_1','slice3_2'});
net.addLayer('slice4', dagnn.Slice(), {'upscore_dsn4'},...
            {'slice4_0','slice4_1','slice4_2','slice4_3'});
net.addLayer('slice5', dagnn.Slice(), {'upscore_dsn5'},...
            {'slice5_0','slice5_1','slice5_2','slice5_3','slice5_4'});

% Concat slices corresponding to the same scale ---------------------------
net.addLayer('concat0', dagnn.Concat(),...
            {'slice2_0','slice3_0','slice4_0','slice5_0'}, {'concat0'});
net.addLayer('concat1', dagnn.Concat(),...
            {'slice2_1','slice3_1','slice4_1','slice5_1'}, {'concat1'});
net.addLayer('concat2', dagnn.Concat(),...
            {'slice3_2','slice4_2','slice5_2'}, {'concat2'});
net.addLayer('concat3', dagnn.Concat(),...
            {'slice4_3','slice5_3'}, {'concat3'});
      
% Fuse scores corresponding to the same scale -----------------------------
net.addLayer('cat0_score', dagnn.Conv('size',[1,1,4,1]), ...
            {'concat0'}, {'concat0_score'}, {'cat0_score_w','cat0_score_b'});
net.addLayer('cat1_score', dagnn.Conv('size',[1,1,4,1]), ...
            {'concat1'}, {'concat1_score'}, {'cat1_score_w','cat1_score_b'});
net.addLayer('cat2_score', dagnn.Conv('size',[1,1,3,1]), ...
            {'concat2'}, {'concat2_score'}, {'cat2_score_w','cat2_score_b'});
net.addLayer('cat3_score', dagnn.Conv('size',[1,1,2,1]), ...
            {'concat3'}, {'concat3_score'}, {'cat3_score_w','cat3_score_b'});
net.addLayer('cat4_score', dagnn.Conv('size',[1,1,1,1]), ...
            {'slice5_4'},{'concat4_score'}, {'cat4_score_w','cat4_score_b'});

% Fuse scores from all scales ---------------------------------------------
net.addLayer('concat_fuse', dagnn.Concat(), ...
            {'concat0_score','concat1_score','concat2_score',...
             'concat3_score','concat4_score'}, {'concat_fuse'});

% Final loss layer --------------------------------------------------------
net.addLayer('loss_fuse',dagnn.Loss(),{'concat_fuse','labels'},{'fuse_loss'});
        