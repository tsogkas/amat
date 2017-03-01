function [lambdaOpt,kappaOpt] = computeOptimalSmoothingParameters()
% Set up parameter grid ---------------------------------------------------
lambdaMin = 0.001;
lambdaMax = 0.1;
lambdaStep= 0.002;
kappaMin  = 1.1;
kappaMax  = 2;
kappaStep = 0.1;
lambdas = lambdaMin:lambdaStep:lambdaMax;
kappas  = kappaMin:kappaStep:kappaMax;
stats(numel(lambdas),numel(kappas)) = struct();

% Read validation images --------------------------------------------------
paths = setPaths();
opts.nThresh= 1;
opts.imPath = fullfile(paths.bsds500im,'val');
opts.gtPath = fullfile(paths.bsds500gt,'val');
imageList   = dir(fullfile(opts.imPath, '*jpg'));

for l=1:numel(lambdas)
    for k=1:numel(kappas)
        stats(l,k) = paramStats(imageList,lambdas(l),kappas(k),opts);
        [stats(l,k).odsP,  stats(l,k).odsR, ...
         stats(l,k).odsF,  stats(l,k).odsT, ...
         stats(l,k).oisP,  stats(l,k).oisR, ...
         stats(l,k).oisF,  stats(l,k).AP] = ...
            computeDatasetStats(stats, opts);
    end
end

% Find optimal combination of parameters
idx = min(cat(1, stats(:).odsF));
[idxl,idxk] = ind2sub([numel(lamdbas),numel(kappas)],idx);
lambdaOpt = lambdas(idxl);
kappaOpt  = kappas(idxk);


% -------------------------------------------------------------------------
function stats = paramStats(imageList,lambda,kappa,opts)
% -------------------------------------------------------------------------
% Initialize stats
opts.nImages = numel(imageList);
cntP = zeros(opts.nImages, opts.nThresh);
cntR = zeros(opts.nImages, opts.nThresh);
sumP = zeros(opts.nImages, opts.nThresh);
sumR = zeros(opts.nImages, opts.nThresh);
scores = zeros(opts.nImages, 4); % optimal P,R,F,T for each image

ticStart = tic;
% parfor (i=1:opts.nImages, opts.parpoolSize)
for i=1:opts.nImages % keep that just for debugging
    % Load image and groundtruth data from disk
    [~,iid,~] = fileparts(imageList(i).name);
    tmp = load(fullfile(opts.gtPath,[iid '.mat' ])); tmp = tmp.groundTruth;
    gt  = false([size(tmp{1}.Boundaries), numel(tmp)]);
    for s=1:numel(tmp), gt(:,:,s) = tmp{s}.Boundaries; end
    img = imread(fullfile(opts.imPath,imageList(i).name));
    % Compute edges on smoothed image
    epb = smoothedEdges(img,lambda,kappa);
    [cntP(i,:), sumP(i,:), cntR(i,:), sumR(i,:),scores(i,:)] = ...
        computeImageStats(epb,gt,opts);
    
    msg = sprintf('Testing params: lambda=%.3f, kappa=%.3f', lambda, kappa);
    progress(msg,i,opts.nImages,ticStart,-1);
end

% Store stats in model struct
stats.cntP = cntP;
stats.sumP = sumP;
stats.cntR = cntR;
stats.sumR = sumR;
stats.scores = scores;

% -------------------------------------------------------------------------
function e = smoothedEdges(img,lambda,kappa)
% -------------------------------------------------------------------------
imgSmoothed = L0Smoothing(img,lambda,kappa);
e = edge(imgSmoothed);

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
        if opts.visualize, plotMatch(1,bmap,gt(:,:,s),match1,match2); end
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
else    % binary symmetry maps
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
