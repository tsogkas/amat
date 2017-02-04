function experimentsMedialAxisExtraction(set)
% Only for BSDS500 dataset at this point
if nargin < 1, set = 'val'; end

paths = setPaths();
MAT = load(['./output/MATBSDS500' upper(set) '.mat']); MAT = MAT.MAT;
savePath = '';
for i=1:numel(MAT)
    img = MAT(i).img;
    seg = MAT(i).seg;
    
    % Compute medial points
    a = amat(img);
    b = mil(img);
    c = deepskeleton(img);
    
    % Compute precision, recall, f-measure
    acc.a = evaluate(a,seg);
    acc.b = evaluate(b,seg);
    acc.c = evaluate(c,seg);
    
    % Save results
    save(savePath, 'acc');
end
