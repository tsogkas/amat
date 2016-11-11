%% Greedy approximation of the weighted set cover problem associated with AMAT
%% Initializations
amat.input          = img;
amat.reconstruction = zeros(H*W,numChannels);
amat.axis           = zeros(H,W,numChannels);
amat.radius         = zeros(H,W);
amat.depth          = zeros(H,W); 
amat.covered        = false(H,W);
% Used in the greedy approximate algorithm. Not sure how we will exploit it
amat.price          = inf(H,W); 
                                

% Easy way to compute the number of NEW pixels that will be covered by each 
% disk if it is added in the solution, taking into account the fact that
% larger disks exceed image boundaries.
numNewPixelsCovered = ones(H,W,numScales);
for r=1:numScales
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), double(filters{r}),'same');
end
numNewPixelsCovered = reshape(numNewPixelsCovered ,H*W,numScales);

%% Test error balancing 
% We must add a scale-based regularization term, to favour larger radii
% even when the errors are 0, in which case dividing by the respective
% radius would not change ordering.
lambda = 1e-1;
reconstructionError = bsxfun(@plus, reconstructionError0, lambda./(1:numScales));
costEffectiveness = reconstructionError./ numNewPixelsCovered;
% Sort costs in ascending order and visualize top disks.
[sortedCosts, indSorted] = sort(costEffectiveness(:),'ascend');
top = 1e2;
[yy,xx,rr] = ind2sub([H,W,numScales], indSorted(1:top));
imshow(amat.input); plotCircles([xx,yy],rr,'w');


%% Debug weighting the appearance and radius errors (try to favour bigger radii)
yc = 24; xc = 40; % test point 
close all;
imshow(amat.input)
plotCircles(repmat([xc,yc], [numScales,1]), (1:numScales)', 'k');
alpha = 0.60;
reconstructionError = bsxfun(@plus, alpha*reconstructionError0, (1-alpha)./(1:numScales));
costEffectiveness = reconstructionError./ numNewPixelsCovered;
[minCost, indMin] = min(costEffectiveness(sub2ind([H,W],yc,xc),:));

%%
pn = [xc yc; xc-1 yc; xc+1 yc; xc yc-1; xc yc+1; xc-1 yc-1; xc+1 yc+1; xc-1 yc+1; xc+1 yc-1];
for i=1:size(pn,1)
    canvas = amat.input; 
    canvas(repmat((x-pn(i,1)).^2 + (y-pn(i,2)).^2 <= indMin^2, [1 1 3])) = 1;
    subplot(3,3,i); imshow(canvas); plotCircles([pn(i,1),pn(i,2)], indMin, 'k');
    title(['Cost: '  num2str(costEffectiveness(sub2ind([H,W],pn(i,2),pn(i,1)),indMin))])
end

for i=1:size(pn,1)
    canvas = amat.input; 
    canvas(repmat((x-pn(i,1)).^2 + (y-pn(i,2)).^2 <= indMin^2, [1 1 3])) = 1;
    subplot(2,9,i); imshow(canvas); plotCircles([pn(i,1),pn(i,2)], indMin, 'k');
    title(['Cost: '  num2str(costEffectiveness(sub2ind([H,W],pn(i,2),pn(i,1)),indMin))])
    canvas = amat.input; 
    canvas(repmat((x-pn(i,1)).^2 + (y-pn(i,2)).^2 <= (indMin+1)^2, [1 1 3])) = 1;
    subplot(2,9,i+9); imshow(canvas); plotCircles([pn(i,1),pn(i,2)], indMin+1, 'k');
    title(['Cost(radius+1): '  num2str(costEffectiveness(sub2ind([H,W],pn(i,2),pn(i,1)),indMin+1))])
end

for i=indMin:numScales
    canvas = amat.input; 
    canvas(repmat((x-xc).^2 + (y-yc).^2 <= i^2, [1 1 3])) = 1;
    subplot(4,4,i-indMin+1); imshow(canvas); plotCircles([xc,yc], i, 'k');
    title(['Cost: '  num2str(costEffectiveness(sub2ind([H,W],xc,yc),i)) ' r=' num2str(i)])
end    

%% Run the greedy algorithm
while ~all(amat.covered(:))
    % Find the most cost-effective set in the current iteration
    [minCost, indMin] = min(costEffectiveness(:));
    % Build set D on the fly
    [yc,xc,r] = ind2sub([H,W,numScales], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= r^2;
    newPixelsCovered = D & ~amat.covered;
    
    % Update AMAT
    amat.price(newPixelsCovered) = minCost; % Set price for new elements
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = repmat(f(yc,xc,:,r), [nnz(newPixelsCovered),1]);
    amat.axis(yc,xc,:) = f(yc,xc,:,r);
    amat.radius(yc,xc) = r;

    % TL;DR: Subtract newly covered pixels from all overlapping disks. 
    % LOND VERSION: We have decided to place a disk centered at (xc,yc)
    % with radius r. This disks covers all points (xd,yd) that lie within
    % radius r. For each point (xd,yd) we must find how many disks it is
    % covered by and subtract 1 from the number of NEW pixels covered by
    % these disks (if they were to be selected as part of the solution).
    % Each (xd,yd) is covered by all disks that have centers at distance <=
    % maxRadius from them.
    [yd,xd] = ind2sub([H,W],find(newPixelsCovered)); % find coordinates of all newly covered pixels inside D
    for i=1:numel(yd)
        distFromCoveredPointSquared = (x-xd(i)).^2 + (y-yd(i)).^2;
        % All the points that lie within max radius distance from the 
        % current point (xd,yd) are centers of (at least) one disk that contains it.
        centersOfCoveringDisks = distFromCoveredPointSquared <= numScales^2;
        % All disks of radius >= minCoveringDistance cover (xd,yd) so we 
        % have to substract 1 from the the number of NEW  pixels that will
        % be covered by each disk if it's added in the solution.
%         minCoveringDistanceSquared = ceil(distFromCoveredPointSquared(centersOfCoveringDisks));
        for rr=1:numScales
            inds = centersOfCoveringDisks & (distFromCoveredPointSquared <= rr^2);
            numNewPixelsCovered(inds, rr) = numNewPixelsCovered(inds, rr) - 1;
        end        
    end
    
    % Update cost effectiveness score
    costEffectiveness = reconstructionError./ numNewPixelsCovered;
    assert(all(isinf(costEffectiveness(sub2ind([H,W],yc,xc), 1:r))))

    % Visualizations
    if 0
        figure(2); clf;
        [rowremove,colremove] = ind2sub([H,W], pixel(remove));
        imshow(amat.input);
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
        subplot(221); imshow(amat.input); plotCircles([xc,yc],r, 'k');
        title('Black: maximum medial disk')
        subplot(222); imshow(amat.covered); title(sprintf('Visited pixels %d/%d',nnz(amat.covered),H*W))
        subplot(223); imshow(amat.axis); title('AMAT axes')
        subplot(224); imshow(amat.radius,[]); title('AMAT radii')
    end
end
amat.reconstruction = reshape(amat.reconstruction,H,W,numChannels);

%% Visualize results
figure(3); clf;
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title('Reconstructed image');