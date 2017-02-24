function branches = groupMedialPoints(mat,scales)
if nargin < 2, scales = 1:40; end
R = numel(scales);
[H,W,C] = size(mat.input);

% Compute individual radius maps and connected components
rad = zeros([H,W,R]);
cc  = cell(1,R);
for r=1:R
    rad(:,:,r) = mat.radius == scales(r);
    cc{r} = bwconncomp(rad(:,:,r));
end

% Grouping scheme
% - TODO: Maybe use a criterion that finds the most appropriate group to
%   merge with, in case of multiple matches.
% - TODO: when I merge two groups, I have to also merge children groups
mask = false(H,W);
maxLabel = 1;

% Scales loop -------------------------------------------------------------
for r=1:R 
%     figure(1); imshow2(sum(rad(:,:,1:r-1),3), rad(:,:,r));
    cc{r}.labels = zeros(1, cc{r}.NumObjects); % zero for non-examined ccs
    margin = ceil(0.5*r);
    % CC loop -------------------------------------------------------------
    for i=1:cc{r}.NumObjects;
        % Create rectangle area around cc for efficiency
        [y,x] = ind2sub([H,W], cc{r}.PixelIdxList{i});
        xmin = max(1,min(x)-margin); xmax = min(W,max(x)+margin);
        ymin = max(1,min(y)-margin); ymax = min(H,max(y)+margin);
        % reset submask
        mask(:) = false;
        mask(cc{r}.PixelIdxList{i}) = true;
        mask(ymin:ymax,xmin:xmax) = bwdist(mask(ymin:ymax,xmin:xmax)) <= r+margin;
        
        % The cc is assigned a new label, unless it has already been merged 
        if cc{r}.labels(i) == 0
            cc{r}.labels(i) = maxLabel;
        end
        
        % Find groups at smaller scales that can potentially be merged 
        mergedLabels = cc{r}.labels(i);
        for rr=(r-1):-1:max(1,r-3)
            for j=1:cc{rr}.NumObjects
                if ~any(mergedLabels == cc{rr}.labels(j)) && any(mask(cc{rr}.PixelIdxList{j}))
                    mergedLabels = [mergedLabels, cc{rr}.labels(j)];
                end
            end
        end    
        
        % Merge labels (use the smallest label as the common label)
        commonLabel = min(mergedLabels);
        for rr=1:r
            cc{rr}.labels(ismember(cc{rr}.labels, mergedLabels)) = commonLabel;
        end
                        
        % Merge ccs at the same scale 
        for j=(i+1):cc{r}.NumObjects
            if any(mask(cc{r}.PixelIdxList{j}))
                cc{r}.labels(j) = cc{r}.labels(i);
            end
        end
        
        % If the component has not been merged, increase maxLabel
        if cc{r}.labels(i) == maxLabel
            maxLabel = maxLabel + 1;
        end
        
    end
end

% Construct label map
branches = zeros(H,W,'uint8');
for r=1:R
    for i=1:cc{r}.NumObjects
        branches(cc{r}.PixelIdxList{i}) = cc{r}.labels(i);
    end
end


    function plotMergedBranches(cc1,cc2)
        tmp = zeros(H,W,'uint8');
        tmp(cc1) = 1;
        tmp(cc2) = 2;
        imagesc(tmp); axis off image;
        rectangle('Position',[xmin,ymin,xmax-xmin+1,ymax-ymin+1],'EdgeColor','y')
    end
end