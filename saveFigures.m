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
img = imread(fullfile(figPath, 'shape.png'));
bnd = imdilate(bwperim(img), strel('disk',5)); bnd(border(bnd,5)) = 0;
bnd = cat(3,bnd, zeros([size(bnd,1) size(bnd,2) 2]));
img = im2double(repmat(img,[1 1 3]));
dark = img/2;
comp = max(dark,bnd);
fig = figure; imshow(comp)
export_fig(fullfile(figPath, 'binary_shape.pdf'),'-transparent',fig);

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
iids = {'3096','54082','85048','295087','42049','101087','86016','145086','302008','253055'};
c = 0.4; zerocolorBranches = [c c c];
parfor i=1:1
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
    [~,idxSeg] = max(max(reshape(ex.seg,[],size(ex.seg,3)),[],1));
    imwrite(imgResized, fullfile(figPath, [iid '_resized.jpg']))
%     fig = figure; imshow(label2rgb(mat.branches, parula(max(mat.branches(:))), zerocolor,'shuffle'));
%     export_fig(fullfile(figPath, [iid '_branches.pdf']), '-transparent', fig)
    fig = figure; imshow(label2rgb(matrefined.branches, parula(max(matrefined.branches(:))), zerocolor,'shuffle'));
    export_fig(fullfile(figPath, [iid '_branches_simplified.pdf']), '-transparent', fig)
    fig = figure; imshow(1-ex.pts(:,:,idxSeg)); 
    export_fig(fullfile(figPath, [iid '_gt_skel.pdf']), '-transparent', fig)
%     ax  = matrefined.axis; ax(repmat(all(matrefined.axis==0,3),[1 1 3])) = 0;
    fig = figure; imshow(matrefined.axis)
    export_fig(fullfile(figPath, [iid '_axes_simplified.pdf']), '-transparent', fig)
end
    

%% Image reconstruction results
% iids = {'3096','54082','85048','295087','42049','101087','86016','145086','302008','253055'};
iids = {'175032','175043','189080','253027','302008','291000','300091','145086','302008','101087','106024','119082','208001'};
for i=2:numel(iids)
    iid = iids{i};
    ex = BMAX500.val(strcmp(iid,{BMAX500.val(:).iid}));
    smoothedResized = imresize(L0Smoothing(ex.img),0.5);
    imgResized = imresize(ex.img, 0.5);
    segResized = imresize(ex.seg, 0.5, 'nearest');
    mat = amat(smoothedResized);
    mat.branches = groupMedialPoints(mat);
    matrefined = refineMAT(mat);
    recmat = reshape(inpaint_nans(double(matrefined.reconstruction)), size(imgResized,1), size(imgResized,2), []);
    recgtseg = seg2reconstruction(imgResized,segResized);
    recgtskel= gtskel2reconstruction(imgResized, imresizeCrisp(ex.pts, 0.5), 0.5*imresizeCrisp(ex.rad, 0.5));
    recmil = spbmil2reconstruction(imgResized);
    idxBest = 1;
    SSIM = ssim(double(recgtskel(:,:,:,1)), im2double(imgResized));
    % Find best segmentation
    for s=2:size(recgtskel,4)
        newssim = ssim(double(recgtskel(:,:,:,s)), im2double(imgResized));
        if newssim > SSIM
            SSIM = newssim; 
            idxBest = s;
        end
    end    
    imwrite(imgResized, fullfile(figPath, [iid '_resized.jpg']))
    imwrite(recmil, fullfile(figPath, [iid '_rec_mil.png']))
    imwrite(recgtseg(:,:,:,idxBest), fullfile(figPath, [iid '_rec_gtseg.png']))
    imwrite(recmat, fullfile(figPath, [iid '_rec_amat.png']))
    imwrite(recgtskel(:,:,:,idxBest), fullfile(figPath, [iid '_rec_gtskel.png']))
end

%% Compute reconstruction and compression metrics 
set = 'val';
% AMAT
tmp = load('/home/tsogkas/code/amat/output/models/amat.mat'); 
modelamat = tmp.model;
mseamat = mean(modelamat.BSDS500.(set).stats.mse);
psnramat = mean(modelamat.BSDS500.(set).stats.psnr);
ssimamat = mean(modelamat.BSDS500.(set).stats.ssim);
compamat = 1./mean(modelamat.BSDS500.(set).stats.compression);

% MIL
tmp = load('/home/tsogkas/code/amat/output/models/model-1000-color-nor-balanced-train.mat'); 
modelmil = tmp.model;
msemil = mean(modelmil.BSDS500.(set).stats.mse);
psnrmil = mean(modelmil.BSDS500.(set).stats.psnr);
ssimmil = mean(modelmil.BSDS500.(set).stats.ssim);
compmil = 1./mean(modelmil.BSDS500.(set).stats.compression);

% GTSET
tmp = load('/home/tsogkas/code/amat/output/models/gtseg.mat'); 
modelgtseg = tmp.model;
msegtseg = mean(modelgtseg.BSDS500.(set).stats.mse);
psnrgtseg = mean(modelgtseg.BSDS500.(set).stats.psnr);
ssimgtseg = mean(modelgtseg.BSDS500.(set).stats.ssim);
compgtseg = 1./mean(modelgtseg.BSDS500.(set).stats.compression);






