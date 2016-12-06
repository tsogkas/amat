%% Setup global parameters and preprocess image
R = 40; % #scales
B = 64; % #bins
errorType = 'se';
encodingType = 'hist';
colorWeights = [];

imgRGB = im2double(imresize(imread('google.jpg'), [128 128], 'nearest')); 
% imgRGB = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/66075.jpg'), [128 128], 'nearest')); 
% imgRGB = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/41004.jpg'), [128 128], 'nearest')); 
[H,W,C] = size(imgRGB);
imgLab = rgb2labNormalized(imgRGB);
imgClustered = clusterImageValues(imgLab, 5); % simplify input

%% Construct filters, calculate perimeneters and disk areas
filters = cell(R,1); for r=1:R, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% Compute individual encodings for color and texture
img = cat(3, binImage(imgLab,B),textonMap(imgRGB, B));
f = imageEncoding(img,filters,'hist',B);


%% Compute decodings g and reconstruction errors at all points and scales
g = f; % proxy decodings (no histogram summarization)

reconstructionError = zeros(H,W,C,R);
maximalityError


%% For testing
p = [83,65];
r1 = 17; 
dr = 3;
B  = 32;
figure(1);
subplot(211); bar(patchEncoding(diskPatch(imgRGB,p,r1),'hist',B)); title(['r=' num2str(r1)])
subplot(212); bar(patchEncoding(diskPatch(imgRGB,p,r1+dr),'hist',B)); title(['r=' num2str(r1+dr)])

%% Test black and white patterns
b = 32;
sz = [100,100];
black = zeros(sz);
white = ones(sz);
half = [black(:,1:sz(2)/2), white(:,1:sz(2)/2)];
checker = white; checker(1:2:end-1,1:2:end) = 0; checker(2:2:end,2:2:end) = 0;
checkerinv = 1-checker;
gray = mean(half(:))*ones(sz);
hblack = patchEncoding(binImage(black(:),B),'hist-normalized',B);
hwhite = patchEncoding(binImage(white(:),B),'hist-normalized',B);
hhalf = patchEncoding(binImage(half(:),B),'hist-normalized',B);
hchecker = patchEncoding(binImage(checker(:),B),'hist-normalized',B);
hcheckerinv = patchEncoding(binImage(checkerinv(:),B),'hist-normalized',B);

thalf = textonMap(half,B);
tchecker = textonMap(checker,B);
htexhalf = patchEncoding(thalf(:),'hist-normalized',B);
htexchecker = patchEncoding(tchecker(:),'hist-normalized',B);



