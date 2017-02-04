function smoothDataset(dataset,overwrite,lambda,kappa)

if nargin < 4, kappa = 2.0; end
if nargin < 3, lambda = 2e-2; end
if nargin < 2, overwrite = false; end
if nargin < 1, dataset = 'bsds'; end


switch dataset
    case 'bsds'
        smoothBSDS500(overwrite,lambda,kappa)
    otherwise, error('Dataset is not supported.')
end


function smoothBSDS500(overwrite,lambda,kappa)
paths = setPaths();
sets = dir(fullfile(paths.bsds500,'images'));
sets(1:2) = []; % remove '.' and '..' dirs['Saving smoothed images for BSDS500 (', sets(s).name, ')...']
for s=1:numel(sets)
    pathOriginal = fullfile(paths.bsds500,'images',sets(s).name);
    pathSmoothed = fullfile(paths.bsds500,...
        sprintf('imagesSmoothed,l=%.2g,k=%.2g',lambda,kappa),sets(s).name);
    mkdir(pathSmoothed)
    images = dir(fullfile(pathOriginal,'*.jpg'));
    msg = sprintf('Saving smoothed images for BSDS500 (%s), l=%.3g, k=%.1g...',...
        sets(s).name, lambda, kappa);
    nImages = numel(images);
    ticStart = tic;
    for i=1:nImages
        if ~exist(fullfile(pathSmoothed,images(i).name),'file') || overwrite
            img = imread(fullfile(pathOriginal,images(i).name));
            smo = L0Smoothing(img,lambda,kappa);
            imwrite(smo,fullfile(pathSmoothed,images(i).name));
        end
        progress(msg,i,nImages,ticStart,10);
    end
end
