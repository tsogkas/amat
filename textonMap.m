function tmap = textonMap(img,numTextons)

if nargin<2, numTextons = 64; end

no = 6; ss = 1; ns = 2; sc = sqrt(2); el = 2;
fname = sprintf('unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,numTextons);
textonData = load(fname); % defines fb,tex,tsim
if size(img,3)==3, tmapim = rgb2gray(img); else tmapim = img; end
tmap = assignTextons(fbRun(textonData.fb,tmapim),textonData.tex);
