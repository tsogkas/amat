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
branches = zeros(H,W,'uint8');
mask = false(H,W);
maxLabel = 1;
for r=1:40
%     figure(1); imshow2(sum(rad(:,:,1:r-1),3), rad(:,:,r));
    cc{r}.labels = zeros(1, cc{r}.NumObjects); % zero for non-examined ccs
    margin = ceil(r);
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
            for rr=(r-1):-1:max(1,r-5)  % First merge ccs of smaller radii 
                for j=1:cc{rr}.NumObjects
                    % merge ccs that are close enough. If cc is matched 
                    % with multiple components, these should have the 
                    % same label anyway.
                    if cc{r}.labels(i) == maxLabel && ...
                       any(mask((cc{rr}.PixelIdxList{j})))
%                         figure(2); plotMergedBranches(cc{r}.PixelIdxList{i},cc{rr}.PixelIdxList{j});
                        cc{r}.labels(i) = cc{rr}.labels(j);
                    end
                end
            end
            % If the component has not been merged, increase maxLabel
            if cc{r}.labels(i) == maxLabel
                maxLabel = maxLabel + 1; 
            end
        end
        
        % Now merge ccs at the same scale 
        for j=i+1:cc{r}.NumObjects
            if cc{r}.labels(i) ~= cc{r}.labels(j) && any(mask(cc{r}.PixelIdxList{j}))
%                 figure(2); plotMergedBranches(cc{r}.PixelIdxList{i},cc{r}.PixelIdxList{j});
                cc{r}.labels(j) = cc{r}.labels(i);
            end
        end
    end
    % Update label map
    for i=1:cc{r}.NumObjects
        branches(cc{r}.PixelIdxList{i}) = cc{r}.labels(i);
    end
end

% Nested utility functions ------------------------------------------------
%     function out = isClose(ccidx)
%         [yy,xx] = ind2sub([H,W], ccidx);
%         xcontained = xx >= xmin & xx <= xmax;
%         ycontained = yy >= ymin & yy <= ymax;
%         out = any(xcontained & ycontained);
%     end

    function plotMergedBranches(cc1,cc2)
        tmp = zeros(H,W,'uint8');
        tmp(cc1) = 1;
        tmp(cc2) = 2;
        imagesc(tmp); axis off image;
        rectangle('Position',[xmin,ymin,xmax-xmin+1,ymax-ymin+1],'EdgeColor','y')
    end
end