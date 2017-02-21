function testReconstruction(models,varargin)

% Default testing options ------------------------------------------------
opts = {'dataset',   'BSDS500',...
        'set',       'val',...   % 'val' or 'test'
        'visualize', false,...
       };                        
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths();
if ischar(opts.set) && strcmp(opts.set, 'val')
    imPath    = fullfile(paths.bsds500im,'val');
    imageList = dir(fullfile(imPath, '*jpg'));
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    imPath    = fullfile(paths.bsds500im,'test');
    imageList = dir(fullfile(imPath, '*jpg'));
elseif isstruct(opts.set)
    disp('Data provided in struct form')
    imageList = opts.set;
    if strcmp(opts.dataset, 'BSDS500')
        if numel(imageList) == 100, opts.set = 'val'; else opts.set = 'test'; end
    end
else
    error('set can be ''val'', ''test'', or a struct containing test data')
end
opts.nImages = numel(imageList);

% Load models and initialize stats ----------------------------------------
if ~iscell(models), models = {models}; end
for m=1:numel(models)
    switch lower(models{m})
        case 'amat'
            models{m} = struct('name',models{m});
        otherwise % load MIL or CNN model
            models{m} = loadModelFromMatFile(models{m},paths);
    end
    models.stats.mse = zeros(opts.nImages, 4);
    models.stats.psnr= zeros(opts.nImages, 4);
    models.stats.ssim= zeros(opts.nImages, 4);
end

% Evaluate models on test images and compute approximate reconstructions --
ticStart = tic;
for i=1:opts.nImages
    if isfield(imageList(i), 'isdir')
        % Load image and groundtruth data from disk
        [~,iid,~] = fileparts(imageList(i).name);
        img = im2double(imread(fullfile(imPath,imageList(i).name)));
    else % Read image and groundtruth from struct
        img = imageList(i).img;
        iid = imageList(i).iid;
    end
    
    clear features 
    for m=1:numel(models)
        switch models{m}.name
            case 'amat'
                rec = reconstructionAMAT(img);
            case 'deepskel'
                rec = reconstructionDeepSkel(models{m},img)
            otherwise % MIL 
                rec = reconstructionMIL(models{m},img);
        end
        [models{m}.stats.cntP(i,:), models{m}.stats.sumP(i,:),...
         models{m}.stats.cntR(i,:), models{m}.stats.sumR(i,:),...
         models{m}.stats.scores(i,:)] = computeImageStats(spb,gt,opts);
    end
    msg = sprintf('Testing on BSDS500 %s set\n', opts.set);
    progress(msg,i,opts.nImages,ticStart,1);
end

% Compute dataset-wide stats
for m=1:numel(models)
    [models{m}.stats.oidP,  models{m}.stats.oidR, ...
     models{m}.stats.oidF,  models{m}.stats.oidT, ...
     models{m}.stats.oisP,  models{m}.stats.oisR, ...
     models{m}.stats.oisF,  models{m}.stats.AP] = ...
        computeDatasetStats(models{m}.stats, opts);
    % Create field with dataset-specific stats
    models{m}.(opts.dataset).(opts.set).stats = stats;
    models{m}.(opts.dataset).(opts.set).opts = opts;
    models{m} = rmfield(models{m},'stats');
    % And store results
    modelPath = fullfile(paths.sbpmil.models, models{m}.name);
    model = models{m}; save(modelPath, 'model')
end

% -------------------------------------------------------------------------
function rec = reconstructionAMAT(img)
% -------------------------------------------------------------------------
mat = amat(img);
rec = mat.reconstruction;

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
function rec = reconstructionMIL(model,img)
% -------------------------------------------------------------------------
spb = spbMIL(img, 'featureSet',model.opts.featureSet,'w',model.w);
spb = spb.thin;

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
