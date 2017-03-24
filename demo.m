%% Demo for the AMAT project
google = imresize(imread('google.jpg'),0.5);
mat = AMAT(google);
figure(1), mat.visualize()
fprintf('Nonzero pixels: %d/%d (%.2f%%)\n', nnz(any(mat.axis,3)), ...
    numel(mat.radius),nnz(any(mat.axis,3))/numel(mat.radius)*100)

%% Interactive demo (this is very slow because of the visualizations)
google = imresize(imread('google.jpg'),0.5); 
figure(2), AMAT(google,'vistop',10); % show 10 highest scoring disks in white color

%% Group medial points into branches (branches are color coded)
mat.group();
figure(3), imagesc(mat.branches); axis off image;

%% Simplify the AMAT
matsimple = mat.simplify();
figure(4), imagesc(matsimple.branches); axis off image;

%% Construct BMAX500 annotations and visualize example annotation
BMAX500 = constructBMAX500(); 
idx = 1;
figure(5);
subplot(121); imshow(BMAX500.val(idx).img); title(['Image ' BMAX500.val(idx).iid])
subplot(122); plotDisks(BMAX500.val(idx).img, BMAX500.val(idx).pts(:,:,1), ...
    BMAX500.val(idx).rad(:,:,1),'sample',0.3)
title('Medial axes (green) and medial disks (red)')

%% Example on training MIL symmetry detector on BMAX500 and testing
milmodel = trainMIL('trainSet',BMAX500.train);
models = testSPB('amat', 'testSet', BMAX500.val);
models = testReconstruction({'amat','gtseg','gtskel'},'testSet', BMAX500.val);

