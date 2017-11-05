% Precompute AMATs for all images in VOC2007

% Train, val, trainval, or test
set = 'test'; 

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
parfor i=1:numel(imageIds)
    img = imread(fullfile(paths.voc2007.images, [imageIds{i}, '.jpg']));
    imgSmoothed = L0Smoothing(img);
    imgSmoothedResized = imresize(imgSmoothed, 0.5, 'bilinear');
    mat = AMAT(imgSmoothedResized);
    parsave(fullfile(savePath, ['amat_' imageIds{i}, '.mat']), mat)
end




