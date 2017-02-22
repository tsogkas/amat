function testReconstruction(models,varargin)

% Default testing options ------------------------------------------------
opts = {'dataset',   'BMAX500',...
        'set',       'val',...   % 'val' or 'test'
        'visualize', false,...
        'parpoolSize', feature('numcores')  % set to 0 to run serially
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
    if strcmp(opts.dataset, 'BMAX500')
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

for m=1:numel(models)
    models{m}.(opts.dataset).(opts.set).stats = stats;
    models{m}.(opts.dataset).(opts.set).opts = opts;
    models{m} = rmfield(models{m},'stats');
    % And store results
    modelPath = fullfile(paths.sbpmil.models, models{m}.name);
    model = models{m}; save(modelPath, 'model')
end

% -------------------------------------------------------------------------
function model = reconstructDataset(model,imageList,opts,paths)
% -------------------------------------------------------------------------
switch lower(model)
    case 'amat'
        model = struct('name',model);
    otherwise % load MIL or CNN model
        model = loadModelFromMatFile(model,paths);
end

% Initialize stats
opts.nImages = numel(imageList);
MSE  = zeros(opts.nImages,1);
PSNR = zeros(opts.nImages,1);
SSIM = zeros(opts.nImages,1);

% Evaluate models on test images and compute approximate reconstructions --
modelName = lower(model.name);
ticStart = tic;
% parfor (i=1:opts.nImages, opts.parpoolSize)
for i=1:opts.nImages
    if isfield(imageList(i), 'isdir')
        img = imread(fullfile(opts.imPath,imageList(i).name));
    else 
        img = imageList(i).img;
    end
    
    switch modelName
        case 'amat'
            rec = reconstructionAMAT(img);
        case 'deepskel'
            rec = reconstructionDeepSkel(model,img);
        otherwise % MIL
            rec = reconstructionMIL(model,img);
    end
    MSE(i)  = immse(rec,im2double(img));
    PSNR(i) = psnr(rec,im2double(img));
    SSIM(i) = ssim(rec, im2double(img));
    msg = sprintf('Testing on %s %s set\n', opts.dataset, opts.set);
    progress(msg,i,opts.nImages,ticStart,-1);
end

models.stats.mse  = MSE;
models.stats.psnr = PSNR;
models.stats.ssim = SSIM;

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
histf = computeHistogramFeatures(img);
spb = spbMIL(img, 'featureSet',model.opts.featureSet,'w',model.w,'histFeatures',histf);
% Get binarized map, scales and thetas
points = spb.thin > model.BMAX500.val.stats.odsT;
scales = spb.scalesMap(points);
thetas = spb.orientMap(points);
% Create rectangular filters at all scales and orientations
filters = cell(numel(histf.thetas), numel(histf.scales));
for s=1:numel(histf.scales)
    sc = histf.scales(s);
    rect = ones(2*sc+1, 2*histf.opts.ratio*sc+1);
    for o=1:numel(histf.thetas)
        filters{o,s} = imrotate(rect,rad2deg(histf.thetas(o)));
    end
end

% Assemble means and create reconstruction
rec = zeros(size(img,1),size(img,2));
[y,x] = find(points);
for i=1:nnz(points)
    yy = y(i); xx = x(i);
    sc = scales(yy,xx);
    th = thetas(yy,xx);
    rf = filters{th,sc};
    [h,w] = size(rf);
    y1 = ceil((h-1)/2); y2 = yy+hs;
    x1 = ceil((w-1)/2); x2 = xx+ws;
    mval = sum(sum(img(y1:y2,x1:x2,:)))/(h*w);
end

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
