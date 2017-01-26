function paths = setPaths()

paths.root = fileparts(mfilename('fullpath')); 

% data
paths.data             = fullfile(paths.root, 'data');
paths.output           = fullfile(paths.root, 'output');
paths.bsds500          = fullfile(paths.data,'BSDS500');
paths.bsds500gtTrain   = fullfile(paths.bsds500,'groundtruth','train');
paths.bsds500gtTest    = fullfile(paths.bsds500,'groundtruth','test');
paths.bsds500gtVal     = fullfile(paths.bsds500,'groundtruth','val');
paths.bsds500imTrain   = fullfile(paths.bsds500,'images','train');
paths.bsds500imTest    = fullfile(paths.bsds500,'images','test');
paths.bsds500imVal     = fullfile(paths.bsds500,'images','val');
