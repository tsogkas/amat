function models = testBoundaryDetection(models,varargin)

% Default testing options ------------------------------------------------
if nargin < 1, models = 'amat'; end
opts = {'dataset',      'BSDS500',...
        'set',          'val',...   % 'val' or 'test'
        'visualize',    false,...
        'parpoolSize',  feature('numcores'),...% set to 0 to run serially
        'edgeDepth',    0,...       % threshold used to get edges from mat
        'cannyThresh',  0,...       % threshold used to get edges from seg
        'minCoverage',  1,...     % threshold used to keep more important segments
        'minSegment',   0,...       % threshold used to keep more important segments
        'nThresh',      30,...      % #thresholds used for computing p-r
        'maxDist',      0.01        % controls max distance of an accurately 
       };                           % detected point from groundtruth.
opts = parseVarargin(opts,varargin,'struct');

% Read test images --------------------------------------------------------
paths = setPaths(); mkdir(paths.amat.models);
if ischar(opts.set) && strcmp(opts.set, 'val')
    opts.imPath = fullfile(paths.bsds500im,'val');
    opts.gtPath = fullfile(paths.bsds500gt,'val');
    imageList   = dir(fullfile(opts.imPath, '*jpg'));
elseif ischar(opts.set) && strcmp(opts.set, 'test')
    opts.imPath = fullfile(paths.bsds500im,'test');
    opts.gtPath = fullfile(paths.bsds500gt,'test');
    imageList   = dir(fullfile(opts.imPath, '*jpg'));
else
    error('set can be ''val'', ''test'', or a struct containing test data')
end

% Load models and initialize stats ----------------------------------------
opts.thresh = linspace(1/(opts.nThresh+1),1-1/(opts.nThresh+1),opts.nThresh)';
if ~iscell(models), models = {models}; end
for m=1:numel(models)
    models{m} = evaluateModel(models{m},imageList,opts,paths);
end

% Compute dataset-wide stats
for m=1:numel(models)
    [models{m}.stats.odsP,  models{m}.stats.odsR, ...
     models{m}.stats.odsF,  models{m}.stats.odsT, ...
     models{m}.stats.oisP,  models{m}.stats.oisR, ...
     models{m}.stats.oisF,  models{m}.stats.AP] = ...
        computeDatasetStats(models{m}.stats, opts);
    % Create field with dataset-specific stats
    models{m}.(opts.dataset).(opts.set).stats = models{m}.stats;
    models{m}.(opts.dataset).(opts.set).opts = opts;
    models{m} = rmfield(models{m},'stats');
    % And store results
    modelPath = fullfile(paths.amat.models, opts2fileName(models{m},opts));
    model = models{m}; save(modelPath, 'model')
end

% -------------------------------------------------------------------------
function model = evaluateModel(model,imageList,opts,paths)
% -------------------------------------------------------------------------
switch lower(model)
    case 'amat'
        opts.thresh = 0.5; opts.nThresh = 1;
        model = loadModelFromMatFile(model,paths);
    otherwise, error('Model not supported')
end

% Initialize stats
opts.nImages = numel(imageList);
cntP = zeros(opts.nImages, opts.nThresh);
cntR = zeros(opts.nImages, opts.nThresh);
sumP = zeros(opts.nImages, opts.nThresh);
sumR = zeros(opts.nImages, opts.nThresh);
scores = zeros(opts.nImages, 4); % optimal P,R,F,T for each image

modelName = lower(model.name);
ticStart = tic;
% parfor (i=1:opts.nImages, opts.parpoolSize)
for i=3:opts.nImages % keep that just for debugging
    % Load image and groundtruth data from disk
    [~,iid,~] = fileparts(imageList(i).name);
    tmp = load(fullfile(opts.gtPath,[iid '.mat' ])); tmp = tmp.groundTruth;
    gt  = false([size(tmp{1}.Boundaries), numel(tmp)]);
    for s=1:numel(tmp), gt(:,:,s) = tmp{s}.Boundaries; end
    
    % Compute edges from 
    epb = amatEdges(['bsds500-' imageList(i).name], opts);
    if size(epb,1) ~= size(gt,1) || size(epb,2) ~= size(gt,2)
        epb = imresize(epb,[size(gt,1),size(gt,2)],'nearest');
    end
    [cntP(i,:), sumP(i,:), cntR(i,:), sumR(i,:),scores(i,:)] = ...
        computeImageStats(epb,gt,opts);
    
    msg = sprintf('Testing %s for boundary detection on %s %s set. ', ...
        modelName, opts.dataset, opts.set);
    progress(msg,i,opts.nImages,ticStart,-1);
end

% Store stats in model struct
model.stats.cntP = cntP;
model.stats.sumP = sumP;
model.stats.cntR = cntR;
model.stats.sumR = sumR;
model.stats.scores = scores;

% -------------------------------------------------------------------------
function epb = amatEdges(img,opts)
% -------------------------------------------------------------------------
mat = amat(img);
if opts.edgeDepth > 0
    epb = mat.depth <= opts.edgeDepth;
elseif opts.cannyThresh > 0
    seg = mat2seg(mat, opts.minCoverage, opts.minCoverage);
    epb = seg2edges(seg, opts.cannyThresh);
else
    seg = mat2seg(mat, opts.minCoverage, opts.minSegment);
    epb = seg2edges(seg)>0;
end

% -------------------------------------------------------------------------
function [cntP,sumP,cntR,sumR,scores] = computeImageStats(pb,gt,opts)
% -------------------------------------------------------------------------
% For Levinstein's method we do not need to threshold
if islogical(pb), 
    thresh = 0.5; pb = double(pb);
else
    thresh = opts.thresh;
end

% Initialize
cntP = zeros(size(thresh));
sumP = zeros(size(thresh));
cntR = zeros(size(thresh));
sumR = zeros(size(thresh));

% Compute numerator (cntX) and denominator (sumX) for precision and recall.
% Note that because there are multiple groundtruth maps to compare with a
% single machine-generated response, the number of true positives, false
% positives etc. used for computing precision is different than the ones
% that are used for computing recall.
for t = 1:numel(thresh),
    % Threshold probability map and thin to 1-pixel width.
    bmap = (pb >= thresh(t));
    bmap = bwmorph(bmap,'thin',Inf);
    
    % Compute matches between symmetry map and all groundtruth maps
    accP = 0;
    for s=1:size(gt,3)
        [match1,match2] = correspondPixels(double(bmap),double(gt(:,:,s)),opts.maxDist);
        if opts.visualize
            plotMatch(1,bmap,gt(:,:,s),match1,match2); drawnow;
        end
        % accumulate machine matches
        accP = accP | match1;
        cntR(t) = cntR(t) + nnz(match2>0); % tp (for recall)
    end
    cntP(t) = nnz(accP); % tp (for precision)
    sumP(t) = nnz(bmap); % tp + fp (for precision)
    sumR(t) = nnz(gt);   % tp + fn (for recall)
end

% Compute precision (P), recall (R) and f-measure (F). 
P = cntP ./ max(eps, sumP);
R = cntR ./ max(eps, sumR);

% Use linear interpolation to find best P,R,F combination.
% scores contains the image-specific optimal P,R,F, after computing optimal
% thresholds using linear interpolation.
[bestP,bestR,bestF,bestT] = findBestPRF(P,R,thresh);
scores = [bestP,bestR,bestF,bestT];

% -------------------------------------------------------------------------
function [odsP, odsR, odsF, odsT, oisP, oisR, oisF, AP] = computeDatasetStats(stats,opts)
% -------------------------------------------------------------------------
% Two standard F-based performance metrics are computed:
% i)  ODS: F-measure for a dataset-wide specified optimal threshold.
% ii) OIS: F-measure for an image-specific optimal threshold.
% iii)AP:  Average precision - equivalent to AUC (area under curve).

P = sum(stats.cntP,1) ./ max(eps, sum(stats.sumP,1));
R = sum(stats.cntR,1) ./ max(eps, sum(stats.sumR,1));

if length(P) > 1 % soft probability maps
    % ODS scores (scalars)
    [odsP,odsR,odsF,odsT] = findBestPRF(P,R,opts.thresh);

    % OIS scores (scalars)
    Pi = stats.cntP ./ max(eps, stats.sumP);
    Ri = stats.cntR ./ max(eps, stats.sumR);
    [~,indMaxF] = max(fmeasure(Pi,Ri),[],2);
    indMaxF = sub2ind(size(Pi), (1:size(Pi,1))', indMaxF);
    oisP = sum(stats.cntP(indMaxF)) ./ max(eps, sum(stats.sumP(indMaxF)));
    oisR = sum(stats.cntR(indMaxF)) ./ max(eps, sum(stats.sumR(indMaxF)));
    oisF = fmeasure(oisP,oisR);

    % AP score (scalar)
    AP = interp1(R,P, 0:0.01:1); 
    AP = sum(AP(~isnan(AP)))/100;
else    % binary maps
    odsP = P; odsR = R; odsF = fmeasure(P,R); odsT = 0.5;
    oisP = odsP; oisR = odsR; oisF = odsF; AP = odsP;
end

% -------------------------------------------------------------------------
function F = fmeasure(P,R), F = 2 .* P .* R ./ max(eps, P+R);

% -------------------------------------------------------------------------
function [bestP,bestR,bestF,bestT] = findBestPRF(P,R,T)
% -------------------------------------------------------------------------
if numel(T) == 1
    bestT = T; bestR = R; bestP = P; bestF = fmeasure(P,R); return
end

bestF = -1;
a = linspace(0,1,100); b = 1-a;
for t = 2:numel(T)
    Rt = a.*R(t) + b.*R(t-1);
    Pt = a.*P(t) + b.*P(t-1);
    Tt = a.*T(t) + b.*T(t-1);
    Ft = fmeasure(Pt,Rt);
    [f,indMax] = max(Ft); 
    if f > bestF
        bestF = f; bestT = Tt(indMax);
        bestP = Pt(indMax); bestR = Rt(indMax); 
    end
end

% -------------------------------------------------------------------------
function model = loadModelFromMatFile(model,paths)
% -------------------------------------------------------------------------
if exist(fullfile(paths.amat.models, model),'file') || ...
   exist([fullfile(paths.amat.models, model) '.mat'],'file')
    tmp = load(fullfile(paths.amat.models, model)); model = tmp.model;
else
    model = struct('name',model); 
end

% -------------------------------------------------------------------------
function fileName = opts2fileName(model,opts)
% -------------------------------------------------------------------------
fileName = model.name;
if opts.edgeDepth > 0
    fileName = [fileName '-edgeDepth=' num2str(opts.edgeDepth)];
else
    if opts.minCoverage < 1
        fileName = [fileName '-minCoverage=' num2str(opts.minCoverage)];
    end
    if opts.minSegment > 0
        fileName = [fileName '-minSegment=' num2str(opts.minSegment)];
    end
    if opts.cannyThresh > 0
        fileName = [fileName '-cannyThresh=' num2str(opts.cannyThresh)];
    end
end
fileName = [fileName '.mat'];

