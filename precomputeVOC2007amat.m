function precomputeVOC2007amat(set)

% Create directory to put data 
paths = setPaths();
savePath = fullfile(paths.amat.output, 'precomputed', 'voc2007');
mkdir(savePath);

% Read image ids
fid = fopen(fullfile(paths.voc2007.sets, 'Main', [set, '.txt'])); 
assert(fid > 0);
imageIds = textscan(fid, '%s'); imageIds = imageIds{1};
fclose(fid);

% Compute AMATs
numImages = numel(imageIds);
parfor i=1:numImages
    matfilePath = fullfile(savePath, ['amat_' imageIds{i}, '.mat']);
    if ~exist(matfilePath, 'file')
        fprintf('Computing VOC2007 AMATs, %s set, %d/%d\n', set, i, numImages)
        img = imread(fullfile(paths.voc2007.images, [imageIds{i}, '.jpg']));
        imgSmoothed = L0Smoothing(img);
        imgSmoothedResized = imresize(imgSmoothed, 0.5, 'bilinear');
        mat = AMAT(imgSmoothedResized);
        parsave(matfilePath, mat)
    end
end

function parsave(matfilePath, mat)
save(matfilePath, 'mat')





