function paths = setPaths()

paths.amat.root = fileparts(mfilename('fullpath')); 

% data
paths.data             = fullfile(paths.amat.root, 'data');
paths.bsds500          = fullfile(paths.data,'BSDS500');
paths.bsds500gt        = fullfile(paths.bsds500, 'groundtruth');
paths.bsds500im        = fullfile(paths.bsds500, 'images');
paths.symmax500        = fullfile(paths.data, 'SYMMAX500');
paths.vgg16            = fullfile(paths.data, 'imagenet-vgg-verydeep-16.mat');
paths.voc2007          = fullfile(paths.data, 'VOCdevkit');

% amat
paths.amat.output      = fullfile(paths.amat.root, 'output');
paths.amat.models      = fullfile(paths.amat.output, 'models');
paths.amat.precomputed = fullfile(paths.amat.output, 'precomputed');
paths.amat.boxes       = fullfile(paths.amat.output, 'boxes');

% spb-mil
paths.spbmil.root      = fullfile(paths.amat.root, 'external', 'spb-mil');
paths.spbmil.output    = fullfile(paths.spbmil.root, 'output');
paths.spbmil.models    = fullfile(paths.spbmil.output, 'models');
paths.spbmil.spectral  = fullfile(paths.spbmil.output, 'spectral');
paths.spbmil.plots     = fullfile(paths.spbmil.output, 'plots');

% edges
paths.edges.root       = fullfile(paths.amat.root,'external','dollar-edges');