function precomputeAMAT(dataset, set, scales, ws, parpoolSize)
if nargin < 1, dataset = 'BSDS500'; end
if nargin < 2, set = 'val'; end
if nargin < 3, scales = 2:41; end
if nargin < 4, ws  = 1e-4; end
if nargin < 5, parpoolSize = feature('numcores'); end

switch lower(dataset)
    case 'voc2007'
        computeVocAmat(set,scales,ws,parpoolSize);
    case 'bsds500'
        computeBsdsAmat(set,scales,ws,parpoolSize);
    otherwise, error('Dataset not supported')
end

% -------------------------------------------------------------------------
function computeVocAmat(set,scales,ws,parpoolSize)
% -------------------------------------------------------------------------
% Setup paths
paths  = setPaths;
imgDir = fullfile(paths.voc2007, 'VOC2007/JPEGImages/');
imgList=fullfile(paths.voc2007, ['VOC2007/ImageSets/Main/' set '.txt']);
addpath(genpath(fullfile(paths.voc2007, 'VOCcode')));

% get list of image ids
if(~exist(imgList,'file')), error('ids file not found'); end
f=fopen(imgList); ids=textscan(f,'%s %*s'); ids=ids{1}; fclose(f);
matDir = fullfile(paths.amat.precomputed,'voc2007', ...
    ['scales=[' num2str(scales) ']-ws=' num2str(ws)]); 
mkdir(matDir);

% Compute AMAT for all images
numImages = numel(ids);
parfor (i=1:numImages, parpoolSize)
% for i=1:numImages
    savePath = fullfile(matDir, ['amat_' ids{i} '.mat']);
    if ~exist(savePath,'file')
        img = imread(fullfile(imgDir, [ids{i} '.jpg']));
        mat = computeAMAT(img,ws,scales);
        saveInParfor(savePath,mat);
    end
end

% -------------------------------------------------------------------------
function computeBsdsAmat(set,scales,ws,parpoolSize)
% -------------------------------------------------------------------------
% Setup paths and get image list
paths = setPaths;
imgDir = fullfile(paths.bsds500im, set);
imgList= dir(fullfile(imgDir,'*.jpg'));
matDir = fullfile(paths.amat.precomputed,'bsds500', ...
    ['scales=[' num2str(scales) ']-ws=' num2str(ws)]); 
mkdir(matDir);

% Compute AMAT for all images
numImages = numel(imgList);
% parfor (i=1:numImages, parpoolSize)
for i=1:numImages
    [~,iid] = fileparts(imgList(i).name);
    savePath = fullfile(matDir, ['amat_' iid '.mat']);
    if ~exist(savePath,'file')
        img = imread(fullfile(imgDir, imgList(i).name));
        mat = computeAMAT(img,scales,ws);
        saveInParfor(savePath,mat);
    end    
end

% -------------------------------------------------------------------------
function mat = computeAMAT(img,scales,ws)
% -------------------------------------------------------------------------
smoothed = L0Smoothing(imresize(img,0.5,'bilinear'));
mat = amat(smoothed,scales,ws);
mat.visualize = []; % remove function handle to save space

function saveInParfor(fileName,mat), save(fileName,'mat');