function img = binImage(img,numBins)
img = max(1, ceil(bsxfun(@times, img, reshape(numBins,1,1,[]))));
