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
VOCxmlPath = fullfile(paths.voc2007.xml,'%s.xml');


%% Extract boxes
parpoolOpen(4);
disp(['Boxes smaller than ' num2str(minBoxWidth) ' pixels will be removed.']);
disp('Obtaining boxes for Pascal 2007 test set...');
numImages = numel(testIms); 
boxesSS = cell(numImages, 1);
boxesAMAT = cell(numImages, 1);

parfor i=1:numImages
    fprintf('Computing boxes for image %d/%d.\n', i, numImages);
    
    % Compute SS boxes
    boxesT    = {};
    priorityT = {};
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

% Filter boxes
for i=1:numel(boxesSS)
    boxesSS{i} = FilterBoxesWidth(boxesSS{i}, minBoxWidth);
    boxesSS{i} = BoxRemoveDuplicates(boxesSS{i});
end
for i=1:numel(boxesAMAT)
    boxesAMAT{i} = FilterBoxesWidth(boxesAMAT{i}, minBoxWidth);
    boxesAMAT{i} = BoxRemoveDuplicates(boxesAMAT{i});
end

save('output/boxes/boxesSS.mat', 'boxesSS');
save('output/boxes/boxesAMAT.mat', 'boxesAMAT')

%%
if ~exist('boxesSS', 'var'), load('output/boxes/boxesSS.mat'); end
if ~exist('boxesAMAT', 'var'), load('output/boxes/boxesAMAT.mat'); end
visualize = 0;
threshIOU = 0.5;

numImages = numel(testIms); 
numBoxesTotal = 0;
numBoxesFoundSS = 0;
numBoxesFoundAMAT = 0;
ticStart = tic;
for i=1:numImages
    rec = VOCreadrecxml(sprintf(VOCxmlPath, testIms{i}));
    img = imread(sprintf(VOCImgPath, testIms{i}));
    bboxGT = cat(1, rec.objects(:).bbox);
    numBoxesGT = size(bboxGT,1);
    numBoxesTotal = numBoxesTotal + numBoxesGT;
    
    % Compute IOU scores of all box proposals and GT boxes.
    iouScoresSS   = zeros(size(boxesSS{i}, 1), numBoxesGT);
    iouScoresAMAT = zeros(size(boxesAMAT{i}, 1), numBoxesGT);
    for b=1:numBoxesGT
        iouScoresSS(:,b)    = iou(boxesSS{i},   bboxGT(b,:));
        iouScoresAMAT(:,b)  = iou(boxesAMAT{i}, bboxGT(b,:));
    end
    
    % Find proposals with max IOU score
    [iouBestSS,  idxBestSS]   = max(iouScoresSS,   [], 1);
    [iouBestAMAT,idxBestAMAT] = max(iouScoresAMAT, [], 1);
    
    numBoxesFoundSS = numBoxesFoundSS + nnz(iouBestSS > 0.5);
    numBoxesFoundAMAT = numBoxesFoundAMAT + nnz(iouBestAMAT > 0.5);
    
    % Visualize results
    if visualize
        imshow(img);
        drawBoxes(bboxGT, 'color','b');
        drawBoxes(boxesSS{i}(idxBestSS,:), 'color', 'm');
        drawBoxes(boxesAMAT{i}(idxBestAMAT,:), 'color', 'y');
    end
    progress('Computing box overlaps...', i, numImages,ticStart, 3);
end


disp('Evaluating the SS boxes on Pascal 2007...')
[boxAboSS, boxMaboSS, boScoresSS, avgNumBoxesSS] = ...
    BoxAverageBestOverlap(gtBoxes, gtImIds, boxesSS);

disp('Evaluating the AMAT boxes on Pascal 2007...')
[boxAboAMAT, boxMaboAMAT, boScoresAMAT, avgNumBoxesAMAT] = ...
    BoxAverageBestOverlap(gtBoxes, gtImIds, boxesAMAT);


fprintf('Mean Average Best Overlap for the SS box-based locations: %.3f\n', boxMaboSS);
fprintf('Mean Average Best Overlap for the AMAT box-based locations: %.3f\n', boxMaboAMAT);
