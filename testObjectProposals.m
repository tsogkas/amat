function method = testObjectProposals(method,varargin)

% Default testing options ------------------------------------------------
if nargin < 1, method = 'edges'; end
opts = {'dataset',      'VOC2007',...
        'set',          'val',...   % 'val' or 'test'
        'visualize',    false,...
        'parpoolSize',  feature('numcores'),...% set to 0 to run serially
        'minCoverage',  0.9         % threshold used to keep more important segments
       };                           
opts = parseVarargin(opts,varargin,'struct');

% Setup paths
paths  = setPaths;
imgDir = fullfile(paths.voc2007, 'VOC2007/JPEGImages/');
imgList=fullfile(paths.voc2007, ['VOC2007/ImageSets/Main/' opts.set '.txt']);
addpath(genpath(fullfile(paths.voc2007, 'VOCcode')));
mkdir(paths.amat.boxes);

% get list of image ids
if(~exist(imgList,'file')), error('ids file not found'); end
f=fopen(imgList); ids=textscan(f,'%s %*s'); ids=ids{1}; fclose(f);

% Load VOC2007 groundtruth boxes
data=boxesData('resDir',[paths.amat.boxes filesep],...
    'dataDir',[paths.voc2007 filesep],'split',opts.set);

% Load structured edges model
model=load(fullfile(paths.edges.root, 'models/forest/modelBsds')); 
model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;

% Load default edge boxes parameters
ebopts = edgeBoxes;
ebopts.alpha = .65;     % step size of sliding window search
ebopts.beta  = .75;     % nms threshold for object proposals
ebopts.minScore = .01;  % min score of boxes to detect
ebopts.maxBoxes = 1e4;  % max number of boxes to detect
ebopts.method = method; 
switch method
    case 'amat-edges'
        nm='AmatEdgeBoxes';
    case 'AmatAxesBoxes';
        nm='AmatAxes';
    case 'edges'
        nm='EdgeBoxes70'; 
    otherwise, error('Method not supported.')
end
ebopts.name=fullfile(paths.amat.boxes, [nm '-' opts.set '.mat']);
edgeBoxes(data.imgs,model,ebopts); ebopts.name=[];
boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);
