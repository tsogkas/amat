%% Modified Selective Search demo on VOC2007 
% Parameters. Note that this controls the number of hierarchical
% segmentations which are combined.
colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};

% Here you specify which similarity functions to use in merging
simFunctionHandles = {@SSSimColourTextureSizeFillOrig, @SSSimTextureSizeFill, @SSSimBoxFillOrig, @SSSimSize};

% Thresholds for the Felzenszwalb and Huttenlocher segmentation algorithm.
% Note that by default, we set minSize = k, and sigma = 0.8.
ks = [50 100 150 300]; % controls size of segments of initial segmentation. 
sigma = 0.8;

% After segmentation, filter out boxes which have a width/height smaller
% than minBoxWidth (default = 20 pixels).
minBoxWidth = 20;

% Comment the following three lines for the 'quality' version
colorTypes = colorTypes(1:2); % 'Fast' uses HSV and Lab
simFunctionHandles = simFunctionHandles(1:2); % Two different merging strategies
ks = ks(1:2);

% Load ground truth boxes and images and image names
load('GroundTruthVOC2007test.mat'); % Loads gtBoxes, gtImIds, gtIms, testIms
paths = setPaths();
VOCImgPath = fullfile(paths.voc2007.images,'%s.jpg');


%% Extract boxes
disp(['Boxes smaller than ' num2str(minBoxWidth) ' pixels will be removed.']);
disp('Obtaining boxes for Pascal 2007 test set...');
numImages = numel(testIms); 
boxesSS = cell(numImages, 1);
boxesAMAT = cell(numImages, 1);

for i=1:numImages
    fprintf('%d ', i);
    
    % Compute SS boxes
    img = imread(sprintf(VOCImgPath, testIms{i}));
    idx = 1;
    for j=1:length(ks)
        k = ks(j); % Segmentation threshold k
        minSize = k; % We set minSize = k
        for n = 1:length(colorTypes)
            colorType = colorTypes{n};
            [boxesT{idx}, ~, ~, ~, priorityT{idx}] = ...
                Image2HierarchicalGrouping(img,sigma,k,minSize,colorType,simFunctionHandles);
            idx = idx + 1;
        end
    end
    
    % Concatenate boxes from all hierarchies
    boxesSS{i} = cat(1, boxesT{:}); 
    
    % Concatenate priorities
    priority = cat(1, priorityT{:}); 
    
    % Do pseudo random sorting as in paper
    priority = priority .* rand(size(priority));
    [~, sortIds] = sort(priority, 'ascend');
    boxesSS{i} = boxesSS{i}(sortIds,:);
    
    % Compute AMAT boxes and do the same
    [blobIndIm, blobBoxes, neighbors] = amat2blobs(testIms{i}, [size(img,1), size(img,2)]);
    [boxesAMAT{i}, ~, priority] = ...
        seg2hierarchicalGrouping(im2double(img), blobIndIm, blobBoxes, neighbors, simFunctionHandles);
    priority = priority .* rand(size(priority));
    [~, sortIds] = sort(priority, 'ascend');
    boxesAMAT{i} = boxesAMAT{i}(sortIds,:);

end
fprintf('\n');

% Filter boxes
for i=1:numel(boxesSS)
    boxesSS{i} = FilterBoxesWidth(boxesSS{i}, minBoxWidth);
    boxesSS{i} = BoxRemoveDuplicates(boxesSS{i});
end
for i=1:numel(boxesAMAT)
    boxesAMAT{i} = FilterBoxesWidth(boxesAMAT{i}, minBoxWidth);
    boxesAMAT{i} = BoxRemoveDuplicates(boxesAMAT{i});
end



%%
disp('Evaluating the SS boxes on Pascal 2007...')
[boxAboSS, boxMaboSS, boScoresSS, avgNumBoxesSS] = ...
    BoxAverageBestOverlap(gtBoxes, gtImIds, boxesSS);

disp('Evaluating the AMAT boxes on Pascal 2007...')
[boxAboAMAT, boxMaboAMAT, boScoresAMAT, avgNumBoxesAMAT] = ...
    BoxAverageBestOverlap(gtBoxes, gtImIds, boxesAMAT);


fprintf('Mean Average Best Overlap for the SS box-based locations: %.3f\n', boxMaboSS);
fprintf('Mean Average Best Overlap for the AMAT box-based locations: %.3f\n', boxMaboAMAT);
