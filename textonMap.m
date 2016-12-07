function tmap = textonMap(img,numTextons,method,numOrient,startSigma,numScales,scaling,elong)

if nargin<2, numTextons = 32; end
if nargin<3, method = 'global'; end
if nargin<4, numOrient = 6; end
if nargin<5, startSigma = 1; end
if nargin<6, numScales = 1; end
if nargin<7, scaling = sqrt(2); end
if nargin<8, elong = 3; end

switch method
    case 'global'
        no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
        fname = sprintf('unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,numTextons);
        textonData = load(fname); % defines fb,tex,tsim
        if size(img,3)==3, tmapim = rgb2gray(img); else tmapim = img; end
        tmap = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);
    case 'image'
        tmap = computeTextons(fbRun(...
            fbCreate(numOrient,startSigma,numScales,scaling,elong),img),numTextons);
    otherwise, error('Method not supported')
end
        
