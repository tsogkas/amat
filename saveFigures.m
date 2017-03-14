% fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, 'test.pdf','-dpdf','-opengl');
figPath = fullfile('report','figures'); mkdir(figPath);
paths = setPaths();
BMAX500 = constructBMAX500();
warning off
%% Teaser
iid = '42049';
ex  = BMAX500.val(strcmp(iid,{BMAX500.val(:).iid}));
mat = amat('bsds500-42049.jpg');
mat.branches = groupMedialPoints(mat);
matrefined = refineMAT(mat);
matrefined.axis(repmat(all(matrefined.axis == 0, 3),[1,1,3])) = 1;
%%
i = 1;
numLabels = numel(unique(ex.seg(:,:,i)));
[H,W,C] = size(ex.img);
fig = figure; imshow(ex.img); 
export_fig(fullfile(figPath, 'teaser_img.pdf'),'-transparent',fig);
fig = figure; plotDisks(label2rgb(ex.seg(:,:,i),...
    parula(numLabels),[0.5 0.5 0.5],'shuffle'),...
    ex.pts(:,:,i),double(ex.rad(:,:,i)).*double(ex.pts(:,:,i)),'sample',0.01);
export_fig(fullfile(figPath, 'teaser_mat.pdf'),'-transparent',fig);
fig = figure; imshow(imresize(matrefined.axis,[H,W],'nearest'));
export_fig(fullfile(figPath, 'teaser_amat.pdf'),'-transparent',fig);
fig  = figure; imshow(imresize(mat.reconstruction,[H,W],'bilinear'));
export_fig(fullfile(figPath, 'teaser_recon.pdf'),'-transparent',fig);

%% Google logo and reconstruction errors
img = imread('google.jpg');


%% Smoothing
iid = '101085'; set = 'val';
imgName = ['bsds500-' iid '.jpg'];
ex  = BMAX500.(set)(strcmp(iid,{BMAX500.(set)(:).iid}));
[H,W,~] = size(ex.img);
smoothed = L0Smoothing(ex.img); 
imwrite(smoothed,fullfile(figPath, [iid '_smoothed.jpg']))

ws = 1e-5; mat = amat(imgName,2:41,ws);
imwrite(imresize(mat.reconstruction,[H,W]),fullfile(figPath, [iid '_recon' num2str(ws) '.jpg']))

ws = 1e-4; mat = amat(imgName,2:41,ws);
imwrite(imresize(mat.reconstruction,[H,W]),fullfile(figPath, [iid '_recon' num2str(ws) '.jpg']))

ws = 1e-3; mat = amat(imgName,2:41,ws);
imwrite(imresize(mat.reconstruction,[H,W]),fullfile(figPath, [iid '_recon' num2str(ws) '.jpg']))

ws = 1e-2; mat = amat(imgName,2:41,ws);
imwrite(imresize(mat.reconstruction,[H,W]),fullfile(figPath, [iid '_recon' num2str(ws) '.jpg']))

ws = 1e-1; mat = amat(imgName,2:41,ws);
imwrite(imresize(mat.reconstruction,[H,W]),fullfile(figPath, [iid '_recon' num2str(ws) '.jpg']))

%% Medial point detection PR-curves
models = {'human','amat','model-1000-color-nor-balanced-train'};
fig = plotPrecisionRecall(models);
export_fig(fullfile(figPath, 'pr.pdf'),'-transparent',fig);

%% Medial point detection qualitative results
paths = setPaths();
iids = {'3096','54082','85048','295087'};
c = 0.4; zerocolor = [c c c];
for i=1:numel(iids)
    iid = iids{i};
    ex  = BMAX500.val(strcmp(iid,{BMAX500.val(:).iid}));
    for s=1:size(ex.seg,3)
        ex.pts(:,:,s) = bwmorph(ex.pts(:,:,s), 'thin',inf);
    end
    ex.pts = imresizeCrisp(ex.pts, 0.5);
    smoothedResized = imresize(L0Smoothing(ex.img),0.5);
    imgResized = imresize(ex.img, 0.5);
    mat = amat(smoothedResized);
    mat.branches = groupMedialPoints(mat);
    matrefined = refineMAT(mat);
    imwrite(imgResized, fullfile(figPath, [iid '_resized.jpg']))
    imwrite(smoothedResized, fullfile(figPath, [iid '_smoothed_resized.jpg']))
    fig = imshow(label2rgb(mat.branches, parula(max(mat.branches(:))), zerocolor,'shuffle'));
    export_fig(fullfile(figPath, [iid '_branches.pdf']), '-transparent', fig)
    fig = imshow(label2rgb(matrefined.branches, parula(max(matrefined.branches(:))), zerocolor,'shuffle'));
    export_fig(fullfile(figPath, [iid '_simplified_branches.pdf']), '-transparent', fig)
    fig = imshow(ex.pts(:,:,1));
end

%% Image reconstruction results
paths = setPaths()
iids = {};
for i=1:numel(iids)
    iids = iids{i};
    ex = BMAX500.val(strcmp(iid,{BMAX500.val(:).iid}));
    smoothedResized = imresize(L0Smoothing(ex.img),0.5);
    imgResized = imresize(ex.img, 0.5);
    mat = amat(smoothedResized);
    mat.branches = groupMedialPoints(mat);
    mat = refineMAT(mat);
    recmat = reshape(inpaint_nans(double(mat.reconstruction)), size(ex.img,1), size(ex.img,2), []);
    
end