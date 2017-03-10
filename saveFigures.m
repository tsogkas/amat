% fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, 'test.pdf','-dpdf','-opengl');
figPath = fullfile('report','figures');
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


