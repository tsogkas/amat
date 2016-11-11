% Greedy search to find medial axis
% TODO: speedup
% TODO: deal with reconstruction from overlapping disks
% TODO: connect medial axis pixels 

% Initialize amat struct
amat.input          = imgClustered;
amat.reconstruction = zeros(H*W,numChannels);
amat.axis           = zeros(H,W,numChannels);
amat.radius         = zeros(H,W);
amat.votes          = zeros(H,W); 
amat.visited        = false(H,W);
% mark the four corners as visited  because they cannot be reached, even by
% a disk with the minimum radius.
amat.visited([1,H,(W-1)*H+1,H*W]) = true; 

tol = 0.05;
scaleSpace      = true(H*W,numScales); % points that are considered as solutions
[rgrid,pgrid]   = meshgrid(1:numScales, 1:H*W);

reconstructionError = reconstructionError0;
% sortedErrors = sort(reconstructionError(:));
while ~all(amat.visited(:))
    % start from the minimum error
%     minError = sortedErrors(find(scaleSpace(:),1));
    minError = min(reconstructionError(:));
    % but also consider points in the scale-space with similar errors
    withineps = abs(reconstructionError - minError) < tol;
    if 0
        figure(1); clf;
        subplot(121); imagesc(reconstructionError); title('Reconstruction error')
        subplot(122); imagesc(withineps); title('Errors within range')
    end
    % and choose the ones corresponding to larger scales first
    [sortedRadii,indSortedRadii] = sort(rgrid(withineps), 'descend');
    [pixel,radius] = find(withineps);
    pixel  = pixel(indSortedRadii); % Kx1 - pairs of (p,r)
    radius = radius(indSortedRadii);
    assert(isequal(radius,sortedRadii))
    while ~isempty(pixel) && ~all(amat.visited(:))
        fprintf('#scale-space candidates in a %.4f-range: %d, #pixels covered: %d/%d\n',...
            tol,numel(pixel),nnz(amat.visited),H*W);
%         disp([' #scale-space candidates within : ' num2str(numel(pixel))]);
        r = radius(1);
        p = pixel(1);

        % mark all pixels that are covered by the disk as visited
        [row,col] = ind2sub([H,W], p);
        distFromCenterSquared = (x-col).^2 + (y-row).^2;
        indisk = distFromCenterSquared <= r^2;
        
        % update AMAT
        amat.votes(indisk) = amat.votes(indisk) + 1;
        amat.visited(indisk) = true;
        amat.reconstruction(indisk,:) = bsxfun(@plus,amat.reconstruction(indisk,:),reshape(f(row,col,:,r),1,3));
%         amat.reconstruction(indisk,:) = repmat(f(row,col,:,r), [nnz(indisk),1]);
        amat.axis(row,col,:) = f(row,col,:,r);
        amat.radius(row,col) = r;
        
        % remove point and all entirely contained disks from the queue
        remove = false(size(pixel));
        for rr=r:-1:1
            fullyContained = ismember(pixel, find(distFromCenterSquared <= (r-rr)^2 & indisk));
            fullyContained = fullyContained & (radius <= rr);
            remove = remove | fullyContained;
            scaleSpace(pixel(fullyContained),1:rr) = false;
            reconstructionError(pixel(fullyContained),1:rr) = inf;
        end
        
        % Visualizations
        if 0
            figure(2); clf;
            [rowremove,colremove] = ind2sub([H,W], pixel(remove));
            imshow(imgClustered);
            plotCircles([colremove,rowremove],radius(remove), 'c'); 
            plotCircles([col,row],r, 'k');
            h1 = plotCircles([col+1 row],r,'m'); title(sprintf('Errors: %f vs %f',reconstructionError0(p,r),reconstructionError0(sub2ind([H,W],row,col+1),r))); 
            delete(h1); h2 = plotCircles([col row+1],r,'m'); title(sprintf('Errors: %f vs %f (diff: %f)',reconstructionError0(p,r),reconstructionError0(sub2ind([H,W],row+1,col),r)));
            delete(h2); h3 = plotCircles([col-1 row],r,'m'); title(sprintf('Errors: %f vs %f (diff: %f)',reconstructionError0(p,r),reconstructionError0(sub2ind([H,W],row,col-1),r)));  
            delete(h3); h4 = plotCircles([col row-1],r,'m'); title(sprintf('Errors: %f vs %f (diff: %f)',reconstructionError0(p,r),reconstructionError0(sub2ind([H,W],row-1,col),r))); 
            delete(h4); plotCircles([col,row],r+1, 'm'); title(sprintf('Errors: %f vs %f (diff: %f)',reconstructionError0(p,r),reconstructionError0(sub2ind([H,W],row,col),r+1))); 
        end
        if 1
            figure(3); clf;
            [rowremove,colremove] = ind2sub([H,W], pixel(remove));
            subplot(221); imshow(imgClustered);
            plotCircles([colremove,rowremove],radius(remove), 'c'); 
            plotCircles([col,row],r, 'k');
            plotCircles([col,row],r+1, 'm');
            title('Black: maximum medial disk | Cyan: fully contained disks')
            subplot(222); imshow(amat.visited); title(sprintf('Visited pixels %d/%d',nnz(amat.visited),H*W))
            subplot(223); imshow(amat.axis); title('AMAT axes')
            subplot(224); imshow(amat.radius,[]); title('AMAT radii')
        end
        assert(remove(1))
        pixel(remove)  = [];
        radius(remove) = [];
    end
end
amat.reconstruction = reshape(amat.reconstruction,H,W,numChannels);
amat.reconstruction = bsxfun(@rdivide, amat.reconstruction, amat.votes);

%% Visualize results
figure(3);
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title('Reconstructed image');
