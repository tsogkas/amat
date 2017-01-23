imageDir = '/home/tsogkas/datasets/BSDS500/images/train';
gtDir = '/home/tsogkas/datasets/BSDS500/groundtruth/train';
imageFiles = dir(fullfile(imageDir,'*.jpg'));
gtFiles = dir(fullfile(gtDir,'*.mat'));

% Plot obtained maximal circles as proof of concept
img = imread(fullfile(imageDir,imageFiles(1).name));
gt  = load(fullfile(gtDir, gtFiles(1).name)); gt = gt.groundTruth;
T = 35; % use a higher threshold to make sure we get good disks
[H,W,C] = size(img);
for s=1:numel(gt)
    seg = gt{s}.Segmentation;
    nSegments = numel(unique(seg));
    pmap = zeros(H,W);
    rmap = zeros(H,W);
    for j=1:nSegments
        segment = seg == j;
        [skel,r] = skeleton(segment);
        p = bwmorph(skel > T,'skel','inf');
        pmap(p) = skel(p);
        rmap(p) = sqrt(r(p));
        if 0
            figure;subplot(121); imshow(pmap,[]); subplot(122); imshow(rmap,[]);
        end
    end
    figure(s); imagesc(seg); axis image off;
    [yy,xx] = find(pmap > 0);
    rr = rmap(pmap > 0);
    assert(numel(xx) == numel(yy) && numel(xx) == numel(rr))
    viscircles([xx,yy],rr,'Color','r','EnhanceVisibility',true,'Linewidth',0.5); 
end
