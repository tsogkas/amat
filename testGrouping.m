%% Compute radius maps and connected components
img = imresize(imread('google.jpg'),0.5,'bilinear');
mat = amat(img);
R = 40;
[H,W,C] = size(mat.input);

rad = zeros([H,W,R]);
for r=1:R
    rad(:,:,r) = mat.radius == r;
end

cc = cell(1,R);
for r=1:R
    cc{r} = bwconncomp(rad(:,:,r));
end

%% 
branches = zeros(H,W,'uint8');
for r=1:R
    for k=1:cc{r}.NumObjects;
        [y,x] = ind2sub([H,W], cc{r}.PixelIdxList{k});
    end
end