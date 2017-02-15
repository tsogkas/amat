function [net, info] = cnnDeepSkeleton(varargin)
% TODO: What should I use as "val" set when I train using trainval?
% TODO: Perhaps weigh the loss of each groundtruth skeleton point based on
%       its skeleton score (confidence).
% TODO: Debug set coverage by maximal disks (should be higher).
% TODO: Add diary/log options.

% Data options
opts.dataDir = constructBMAX500('resize',0.5);
opts.visualizeDataset = 0;
opts.averageImage = [116.66877; 122.67892; 104.00699];  
opts.debug = 1; % use small subset of the data for debugging
opts.mode = 'train';

% Suffix and results directory setup
opts.modelType = 'vgg-vd-16' ;
paths = setPaths();
sfx = opts.modelType ;
opts.expDir = fullfile(paths.amat.output, ['deepskel-' sfx]) ;
[opts, varargin] = vl_argparse(opts, varargin) ;

% Train and GPU options
% Example on how to selecte learning rate:
% In the original DeepSkeleton work they use 10K-20K iterations with Caffe.
% Since in MatConvNet each epoch goes through the complete train set, we
% have to adjust the #epochs to get comparable results.
% - BSDS trainval set: 300 (200 for train set)
% - Augmentation factor: 36 (~36*6 = 216 if we consider each segmentation 
%   as a separate training example).
% - Batchsize: 10
% Hence: #iters/epoch = 300*36 / 10 = 1080, and we need ~15-20 epochs to
% reach 20K iterations (~3 epochs in the case we use all segmentations).
opts.numFetchThreads = 12 ;
opts.train.gpus = 1; 
opts.train.learningRate = [1e-6, 1e-7, 1e-8]; % #epochs = numel(learningRate)
opts.train.momentum = 0.9 ;
opts.train.weightDecay = 0.0002 ;
opts.train.batchSize = 10 ;
[~,machine] = system('hostname'); 
if strcmp(machine,'izual'), opts.train.gpus = []; end
opts = vl_argparse(opts, varargin) ;

% Initialize model
disp('Initializing model...')
rng(0); % set rng for reproducibility
net = cnnInit('initNetPath', paths.vgg16,...
              'learningRate',opts.train.learningRate,...
              'momentum', opts.train.momentum,...
              'weightDecay', opts.train.weightDecay,...
              'batchSize',opts.train.batchSize);

% Prepare training data
disp('Preparing training data...')
imdb = getBMAX500Imdb('dataDir', opts.dataDir,...
                      'mode',opts.mode,...
                      'averageImage', opts.averageImage,...
                      'debug', opts.debug,...
                      'visualizeDataset', opts.visualizeDataset);

% Train
disp('-------------------------------------------------------------------')
disp('Training parameters')
disp('-------------------------------------------------------------------')
[net, info] = cnn_train_dag(net, imdb, @(x,y) getBatch(x,y,opts.train.gpus), ...
                      'expDir', opts.expDir, ...
                      net.meta.trainOpts, ...
                      opts.train) ;


% -------------------------------------------------------------------------
function out = getBatch(imdb, batch, useGpu)
% -------------------------------------------------------------------------
images = single(imdb.images.data(:,:,:,batch));
labels = single(imdb.images.labels(:,:,:,batch));
images = bsxfun(@minus, images, reshape(imdb.meta.averageImage,1,1,3));
if ~isempty(useGpu)
    images = gpuArray(images);
    labels = gpuArray(labels);
end
out = {'input', images, 'label', labels};

% -------------------------------------------------------------------------
function imdb = getBMAX500Imdb(varargin)
% -------------------------------------------------------------------------
opts.mode = 'train';
opts.debug = false;
opts.dataDir = [] ;
opts.lite = false ;
opts.averageImage = zeros(3,1);
opts.visualizeDataset = false;
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
% NOTE: Using the first image to compute Dmax is enough because all images 
% are either DxK or KxD.
Dmax = max(size(bmax500.train(1).img)); 
Hmax = ceil(max(scales)*Dmax); Wmax = Hmax;

% Receptive fields at the level of the scale-associated outputs. Used to
% construct a scale-associated groundtruth map. This is essentially the
% number of labels #labels/L for the equivalent multi-classification
% regression problem
rfields = [14,40,92,196];

% Load only the necessary subsets to save RAM
if strcmp(opts.mode, 'test')
    sets = {'test'};
else
    sets = {'train','val'}; 
end

% Load images/groundtruth and augment scales. Each segmentation is
% considered as an individual training example. 
for s=1:numel(sets)
    set = bmax500.(sets{s});
    % If debug mode, use only one image from each set
    if opts.debug, numImages = 1; else numImages = numel(set); end
    % Total number of groundtruth maps for preallocation
    numGroundtruthTotal = sum(cellfun(@(x) size(x,3), {set(1:numImages).pts}));
    % Preallocate
    imgs = zeros(Hmax,Wmax,3,numGroundtruthTotal*numel(scales),'uint8');
    lbls = zeros(Hmax,Wmax,1,numGroundtruthTotal*numel(scales),'uint8');
    radi = zeros(Hmax,Wmax,1,numGroundtruthTotal*numel(scales),'uint8');
    % Now assemble all examples
    ind = 1;
    for i=1:numImages
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
%             for k=1:size(pts,3)
%                 plotDisks(img,pts(:,:,k),rad(:,:,k));
%             end
            % Center data in the maxH x maxW array
            [H,W,~] = size(img);
            hs = ceil((Hmax-H+1)/2); ws = ceil((Wmax-W+1)/2);
            imgs(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = repmat(img,[1,1,1,ngt]);
            radi(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = reshape(rad,H,W,1,[]);
            lbls(hs:hs+H-1, ws:ws+W-1,:, ind:ind+ngt-1) = reshape(pts,H,W,1,[]);
            ind = ind+ngt;
        end
    end
    % Create scale-associated labels.
    radi = 1.2 * radi; % make sure that the receptive field is large enough
    % MatConvNet ignores 0-labels during the loss computation. We use "1"
    % as the background label.
    lbls(radi == 0) = 1;
    lbls(radi > 0 & radi <= rfields(1)) = 2;
    for r=2:numel(rfields)
        lbls(radi>rfields(r-1) & radi <= rfields(r)) = r+1;
    end
    imdb.(sets{s}).images = imgs;
    imdb.(sets{s}).labels = lbls;
%     imdb.(sets{s}).radius = radi;
end

% Tidy imdb
switch opts.mode
    case 'test'
        imdb.images.data   = imdb.test.images;  imdb.test.images = [];
        imdb.images.labels = imdb.test.labels;  imdb.test.labels = [];
%         imdb.images.radius = imdb.test.radius;  imdb.test.radius = [];
        imdb = rmfield(imdb, 'test');
    case 'train'
        imdb.images.data   = imdb.train.images; imdb.train.images = [];
        imdb.images.labels = imdb.train.labels; imdb.train.labels = []; 
%         imdb.images.radius = imdb.train.radius; imdb.train.radius = []; 
        imdb = rmfield(imdb, 'train');
    case 'trainval'
        imdb.images.data = cat(4, imdb.train.images, imdb.val.images);
        imdb.train.images = []; imdb.val.images = [];
        imdb.images.labels = cat(4, imdb.train.labels, imdb.val.labels);
        imdb.train.labels = []; imdb.val.labels = [];
%         imdb.images.radius = cat(4, imdb.train.radius, imdb.val.radius);
%         imdb.train.radius = []; imdb.val.radius = [];
        imdb = rmfield(imdb, {'train','val'});
    otherwise, error('opts.mode can be ''train'',''trainval'', or''test''')
end

% Visual Inspection of dataset;
if opts.visualizeDataset
    figure(1); montage(imdb.images.data); title('Images')
    figure(2); montage(imdb.images.labels,[1,numel(rfields)+1]); title('Skeletons')
%     figure(3); montage(imdb.images.radius,'DisplayRange',[]); title('Radius maps')
end

% Augment orientations
imgs = imdb.images.data;
lbls = imdb.images.labels;
% radi = imdb.images.radius;
for theta = thetas
    if theta ~= 0
        imdb.images.data   = cat(4, imdb.images.data,   imrotate(imgs,theta,'bilinear'));
        imdb.images.labels = cat(4, imdb.images.labels, imrotate(lbls,theta,'nearest'));
%         imdb.images.radius = cat(4, imdb.images.radius, imrotate(radi,theta,'nearest'));
    end
end
clear imgs lbls radi img pts rad

% Augment flipping
imdb.images.data = cat(4, imdb.images.data,...
    fliplr(imdb.images.data), flipud(imdb.images.data));
imdb.images.labels = cat(4, imdb.images.labels,...
    fliplr(imdb.images.labels), flipud(imdb.images.labels));
% imdb.images.radius = cat(4, imdb.images.radius,...
%     fliplr(imdb.images.radius), flipud(imdb.images.radius));

% Assign set indexes
if strcmp(opts.mode, 'train')
    nTrain = size(imdb.images.data,4);
    nVal = size(imdb.val.images,4);
    imdb.images.data  = cat(4, imdb.images.data,   imdb.val.images);
    imdb.images.labels= cat(4, imdb.images.labels, imdb.val.labels);
    imdb = rmfield(imdb,'val');
    imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal)];
elseif strcmp(opts.mode, 'trainval') 
    nTrain = size(imdb.images.data,4);
    imdb.images.set = [ones(1,nTrain) 2*ones(1,nVal)]; % what is val set ??
else error('Invalid opts.mode')
end
imdb.meta.sets = sets;
imdb.meta.averageImage = opts.averageImage;
imdb.meta.numLabels = numel(rfields) + 1;

% Data sanity checks
dataSize  = size(imdb.images.data);
labelSize = size(imdb.images.labels);
assert(isequal(dataSize([1,2,4]), labelSize([1 2 4])));
assert(all(isinrange(imdb.images.labels,[0,imdb.meta.numLabels])))
assert(all(imdb.images.set <= 3 & imdb.images.set >= 1))

% Display dataset information
disp('-------------------------------------------------------------------')
disp('Training data information')
disp('-------------------------------------------------------------------')
disp(['#images: ' num2str(nTrain) ' (' num2str(dataSize(1:2)) ')'])
trainlbls = imdb.images.labels(:,:,:,imdb.images.set==1);
numPixels = numel(trainlbls);
fprintf('Label stats: %.2f(1), %.2f(2), %.2f(3), %.2f(4), %.2f(5)',...
    nnz(trainlbls==1)/numPixels, nnz(trainlbls==2)/numPixels,...
    nnz(trainlbls==3)/numPixels, nnz(trainlbls==4)/numPixels,...
    nnz(trainlbls==5)/numPixels)

% -------------------------------------------------------------------------
function net = cnnInit(varargin)
% -------------------------------------------------------------------------
opts.initNetPath = [];
opts.cudnnWorkspaceLimit = 1024*1024*1204 ; % 1GB
opts.learningRate = [1e-6, 1e-7, 1e-8] ;
opts.momentum = 0.9 ;
opts.weightDecay = 0.0002 ;
opts.batchSize = 10 ;
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
    net.params(indw).learningRate = 0.05; % 0.01 for cat2-score in deepskel prototxt - probably a typo
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
    net.params(indw).value = val * ones(size(net.params(indw).value),'single');
end

% Keep these layers inputs as they will be used for side outputs
for l={'conv2_2','conv3_3','conv4_3','conv5_3'}
    net.vars(net.getVarIndex(l)).precious = true;
end

% Meta parameters 
net.meta.trainOpts.learningRate = opts.learningRate ;
net.meta.trainOpts.numEpochs = numel(opts.learningRate) ;
net.meta.trainOpts.momentum = opts.momentum ;
net.meta.trainOpts.weightDecay = opts.weightDecay ;
net.meta.trainOpts.batchSize = opts.batchSize ;

% Set objectives that will be considered during backpropagation
net.meta.trainOpts.derOutputs = ...
    {'loss2',1,'loss3',1,'loss4',1,'loss5',1,'loss_fuse',1};

% --------------------------------------------------------------------
function net = vgg16deepskel()
% --------------------------------------------------------------------
% Original HED architecture (https://arxiv.org/abs/1504.06375), extended
% with deep side, loss, and fusion layers for combining skeleton responses
% at different scales.

net = dagnn.DagNN();

% ------------------------------ CONV -------------------------------------
% NOTE: In DeepSkeleton the conv1_1 layer has pad=35. Since we build the
% training images so that most of them have considerable "empty" space
% around the border, we use pad=1. We might need to change that.
% Conv 1 
net.addLayer('conv1_1', dagnn.Conv('size',[3,3,3,64],'stride',1,'pad',35),...
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
         
% --------------------- DSN (Side output layers) --------------------------
% DSN 2
net.addLayer('score_dsn2', dagnn.Conv('size',[1,1,128,2]),...
            {'conv2_2'}, {'score_dsn2'},{'score_dsn2_w', 'score_dsn2_b'});
net.addLayer('deconv2', dagnn.ConvTranspose('size',[4,4,2,2],'upsample',2,'crop',35),...
            {'score_dsn2'}, {'score_dsn2_up'}, {'deconv2_w','deconv2_b'});
% net.addLayer('loss2', dagnn.Loss(),{'score_dsn2_up','label'}, {'loss2'});
net.addLayer('loss2', dagnn.ScaleLoss('scale',2),{'score_dsn2_up','label'}, {'loss2'});
                  
% DSN 3
net.addLayer('score_dsn3', dagnn.Conv('size',[1,1,256,3]),...
            {'conv3_3'}, {'score_dsn3'},{'score_dsn3_w', 'score_dsn3_b'});
net.addLayer('deconv3', dagnn.ConvTranspose('size',[8,8,3,3],'upsample',4,'crop',35),...
            {'score_dsn3'}, {'score_dsn3_up'}, {'deconv3_w','deconv3_b'});
net.addLayer('loss3', dagnn.ScaleLoss('scale',3),{'score_dsn3_up','label'}, {'loss3'});
                  
% DSN 4
net.addLayer('score_dsn4', dagnn.Conv('size',[1,1,512,4]),...
            {'conv4_3'}, {'score_dsn4'},{'score_dsn4_w', 'score_dsn4_b'});
net.addLayer('deconv4', dagnn.ConvTranspose('size',[16,16,4,4],'upsample',8,'crop',35),...
            {'score_dsn4'}, {'score_dsn4_up'}, {'deconv4_w','deconv4_b'});
net.addLayer('loss4', dagnn.ScaleLoss('scale',4),{'score_dsn4_up','label'}, {'loss4'});
         
% DSN 5
net.addLayer('score_dsn5', dagnn.Conv('size',[1,1,512,5]),...
            {'conv5_3'}, {'score_dsn5'},{'score_dsn5_w', 'score_dsn5_b'});
net.addLayer('deconv5', dagnn.ConvTranspose('size',[32,32,5,5],'upsample',16,'crop',39),...
            {'score_dsn5'}, {'score_dsn5_up'}, {'deconv5_w','deconv5_b'});
net.addLayer('loss5', dagnn.ScaleLoss('scale',5),{'score_dsn5_up','label'}, {'loss5'});

% ----------------------- SLICE side outputs ------------------------------
net.addLayer('slice2', dagnn.Slice(), {'score_dsn2_up'},...
            {'slice2_0','slice2_1'}); % scale 0 --> background 
net.addLayer('slice3', dagnn.Slice(), {'score_dsn3_up'},...
            {'slice3_0','slice3_1','slice3_2'});
net.addLayer('slice4', dagnn.Slice(), {'score_dsn4_up'},...
            {'slice4_0','slice4_1','slice4_2','slice4_3'});
net.addLayer('slice5', dagnn.Slice(), {'score_dsn5_up'},...
            {'slice5_0','slice5_1','slice5_2','slice5_3','slice5_4'});

% --------------- CONCAT slices of the same scale -------------------------
net.addLayer('concat0', dagnn.Concat(),...  % scale 0 (background)
            {'slice2_0','slice3_0','slice4_0','slice5_0'}, {'concat0'});
net.addLayer('concat1', dagnn.Concat(),...  % scale 1 
            {'slice2_1','slice3_1','slice4_1','slice5_1'}, {'concat1'});
net.addLayer('concat2', dagnn.Concat(),...  % scale 2
            {'slice3_2','slice4_2','slice5_2'}, {'concat2'});
net.addLayer('concat3', dagnn.Concat(),...  % scale 3 (scale 4 does not need concat)
            {'slice4_3','slice5_3'}, {'concat3'});        
      
% Combine scale-responses from different stages into a single map ---------
net.addLayer('cat0_score', dagnn.Conv('size',[1,1,4,1]), ...    % scale 0 (background)
            {'concat0'}, {'concat0_score'}, {'cat0_score_w','cat0_score_b'});
net.addLayer('cat1_score', dagnn.Conv('size',[1,1,4,1]), ...    % scale 1
            {'concat1'}, {'concat1_score'}, {'cat1_score_w','cat1_score_b'});
net.addLayer('cat2_score', dagnn.Conv('size',[1,1,3,1]), ...    % scale 2
            {'concat2'}, {'concat2_score'}, {'cat2_score_w','cat2_score_b'});
net.addLayer('cat3_score', dagnn.Conv('size',[1,1,2,1]), ...    % scale 3
            {'concat3'}, {'concat3_score'}, {'cat3_score_w','cat3_score_b'});
net.addLayer('cat4_score', dagnn.Conv('size',[1,1,1,1]), ...    % scale 4
            {'slice5_4'},{'concat4_score'}, {'cat4_score_w','cat4_score_b'});

% --------------------- CONCAT scale responses ----------------------------
net.addLayer('concat_fuse', dagnn.Concat(), ...
            {'concat0_score','concat1_score','concat2_score',...
             'concat3_score','concat4_score'}, {'concat_fuse'});

% --------------------- LOSS on fused scores ------------------------------
net.addLayer('loss_fuse',dagnn.ScaleLoss('scale',5),{'concat_fuse','label'},{'loss_fuse'});
        