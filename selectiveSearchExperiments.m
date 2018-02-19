%% Modified Selective Search demo on VOC2007 

% Load ground truth boxes and images and image names
load('GroundTruthVOC2007test.mat'); % Loads gtBoxes, gtImIds, gtIms, testIms
paths = setPaths();
VOCImgPath = fullfile(paths.voc2007.images,'%s.jpg');
VOCxmlPath = fullfile(paths.voc2007.xml,'%s.xml');
fastMode = true;


%% Extract boxes
parpoolOpen(4);
disp('Obtaining boxes for Pascal 2007 test set...');
numImages = numel(testIms); 
boxesSS   = cell(numImages, 1);
boxesAMAT = cell(numImages, 1);

parfor i=1:numImages
    fprintf('Computing boxes for image %d/%d.\n', i, numImages);
    img = imread(sprintf(VOCImgPath, testIms{i}));
    
    boxesSS{i}   = selectiveSearchBoxes(img, fastMode);
    boxesAMAT{i} = amatBoxes(img, testIms{i});

end

save('output/boxes/boxesSS.mat', 'boxesSS');
save('output/boxes/boxesAMAT.mat', 'boxesAMAT')

%%
if ~exist('boxesSS', 'var'), load('output/boxes/boxesSS.mat'); end
if ~exist('boxesAMAT', 'var'), load('output/boxes/boxesAMAT.mat'); end
visualize = 0;
threshIOU = 0.5;

% REMINDER: Selective search uses a [ymin,xmin,ymax,xmax] box format.
numImages = numel(testIms); 
numBoxesTotal = 0;
numBoxesFoundSS = 0;
numBoxesFoundAMAT = 0;
numBoxesFoundCombined = 0;
ticStart = tic;
for i=1:numImages
    rec = VOCreadrecxml(sprintf(VOCxmlPath, testIms{i}));
    bboxGT = cat(1, rec.objects(:).bbox);
    numBoxesGT = size(bboxGT,1);
    numBoxesTotal = numBoxesTotal + numBoxesGT;
    
    % Find best overlaps for each ground-truth box
    [iouBestSS, idxBestSS]     = BoxBestOverlap(bboxGT, boxesSS{i}(:,[2,1,4,3]));
    [iouBestAMAT, idxBestAMAT] = BoxBestOverlap(bboxGT, boxesAMAT{i}(:,[2,1,4,3]));
    
    % Find how many exceed the set IOU threshold
    coveredBySS   = iouBestSS > threshIOU;
    coveredByAMAT = iouBestAMAT > threshIOU;
    numBoxesFoundSS   = numBoxesFoundSS   + nnz(coveredBySS);
    numBoxesFoundAMAT = numBoxesFoundAMAT + nnz(coveredByAMAT);
    numBoxesFoundCombined = numBoxesFoundCombined + nnz(coveredBySS | coveredByAMAT);
    
    % Visualize results
    if visualize
        img = imread(sprintf(VOCImgPath, testIms{i}));
        imshow(img);
%         drawBoxes(bboxGT, 'color','b');
%         drawBoxes(boxesSS{i}(idxBestSS,[2,1,4,3]), 'color', 'm');
%         drawBoxes(boxesAMAT{i}(idxBestAMAT,[2,1,4,3]), 'color', 'y');
        if nnz(xor(coveredByAMAT,coveredBySS))
            drawBoxes(bboxGT((coveredBySS & ~coveredByAMAT) | (~coveredBySS & coveredByAMAT),:), 'color','b');
            drawBoxes(boxesAMAT{i}(idxBestAMAT(coveredBySS & ~coveredByAMAT) ,[2,1,4,3]), 'color', 'r');
            drawBoxes(boxesAMAT{i}(idxBestAMAT(~coveredBySS & coveredByAMAT) ,[2,1,4,3]), 'color', 'g');
            keyboard;
        end
    end
    progress('Computing box overlaps...', i, numImages,ticStart, 3);
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
