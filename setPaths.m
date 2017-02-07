function paths = setPaths()

paths.root = fileparts(mfilename('fullpath')); 

% data
paths.data             = fullfile(paths.root, 'data');
paths.output           = fullfile(paths.root, 'output');
paths.bsds500          = fullfile(paths.data,'BSDS500');
paths.bsds500gt        = fullfile(paths.bsds500, 'groundtruth');
paths.bsds500im        = fullfile(paths.bsds500, 'images');

% TODO: maybe remove the following to reduce clutter?
paths.bsds500gtTrain   = fullfile(paths.bsds500gt,'train');
paths.bsds500gtTest    = fullfile(paths.bsds500gt,'test');
paths.bsds500gtVal     = fullfile(paths.bsds500gt,'val');
paths.bsds500imTrain   = fullfile(paths.bsds500,'images','train');
paths.bsds500imTest    = fullfile(paths.bsds500,'images','test');
paths.bsds500imVal     = fullfile(paths.bsds500,'images','val');
