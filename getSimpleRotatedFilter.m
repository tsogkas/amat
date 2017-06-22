function simpleRotatedFilter =  getSimpleRotatedFilter(filter)
[H,W] = size(filter);
[y,x] = find(filter);
onesInd = [y,x];
simpleRotatedFilter = zeros(H+1,W);
for w=1:W
    wInd = onesInd(onesInd(:,2) == w);
    wMin = min(wInd);
    wMax = max(wInd);
    simpleRotatedFilter(wMin,w) = 1;
    simpleRotatedFilter(wMax+1,w) = -1;
end
end
