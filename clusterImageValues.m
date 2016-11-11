function img = clusterImageValues(img,numClusters)
% Perform k-means to simplify  image

[H,W,numChannels] = size(img);
img = reshape(img, H*W,numChannels);
[clusterInd, clusterCenters] = kmeans(img,numClusters, 'replicates',10);
for k=1:numClusters
    img(clusterInd == k,:) = repmat(clusterCenters(k,:), [nnz(clusterInd == k),1]);
end
img = reshape(img,H,W,numChannels);
