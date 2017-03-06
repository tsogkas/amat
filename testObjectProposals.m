function models = testObjectProposals(models,varargin)

% Default testing options ------------------------------------------------
if nargin < 1, models = 'amat'; end
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
imgList=fullfile(paths.voc2007, ['VOC2007/ImageSets/Main/' set '.txt']);
addpath(genpath(fullfile(paths.voc2007, 'VOCcode')));

% get list of image ids
if(~exist(imgList,'file')), error('ids file not found'); end
f=fopen(imgList); ids=textscan(f,'%s %*s'); ids=ids{1}; fclose(f);

data=boxesData('resDir',paths.amat.boxes,'dataDir',paths.voc2007','split',set);
nm='EdgeBoxes70'; opts.name=fullfile(paths.amat.boxes, [nm '-' set '.mat']);
edgeBoxes(data.imgs,model,opts); opts.name=[];
boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);

% -------------------------------------------------------------------------
function model = loadModelFromMatFile(model,paths)
% -------------------------------------------------------------------------
if exist(fullfile(paths.amat.models, model),'file') || ...
   exist([fullfile(paths.amat.models, model) '.mat'],'file')
    tmp = load(fullfile(paths.amat.models, model)); model = tmp.model;
else
    model = struct('name',model); 
end
