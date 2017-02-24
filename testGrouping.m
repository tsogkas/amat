%% Compute radius maps and connected components
% img = imresize(imread('google.jpg'),0.5,'bilinear');
% mat = amat(img);
branches = groupMedialPoints(mat);
