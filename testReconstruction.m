function testReconstruction()

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
    gtPath    = fullfile(paths.symmax500,'val');
    imageList = dir(fullfile(imPath, '*jpg'));
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    imPath    = fullfile(paths.bsds500im,'test');
    gtPath    = fullfile(paths.symmax500,'test');
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
        case 'levinstein'
            models{m} = loadLevinsteinModel(models{m},paths);
        case 'lindeberg'
            models{m} = struct('name',models{m});
        otherwise % load MIL or CNN model
            models{m} = loadModelFromMatFile(models{m},paths);
    end
end

