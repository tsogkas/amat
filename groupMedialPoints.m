function branches = groupMedialPoints(mat, marginFactor)

if nargin < 2, marginFactor = 0.5; end

% Compute individual radius maps and connected components
R = numel(mat.scales);
for r=R:-1:1
    cc(r) = bwconncomp(mat.radius == mat.scales(r));
end

% Initialize mask and maxLabel
[H,W,C] = size(mat.input);
mask = false(H,W);  % proximity mask
maxLabel = 1;

% For all scales 
for r=1:R 
    cc(r).labels = zeros(1, cc(r).NumObjects); % zero for non-examined ccs
    margin = ceil(marginFactor*r);
    % For all connected components at the same scale
    for i=1:cc(r).NumObjects;
        % Create proximity mask in rectangle around cc for efficiency
        mask(:) = false; mask(cc(r).PixelIdxList{i}) = true;
        [y,x] = ind2sub([H,W], cc(r).PixelIdxList{i});
        xmin = max(1,min(x)-margin); xmax = min(W,max(x)+margin);
        ymin = max(1,min(y)-margin); ymax = min(H,max(y)+margin);
        mask(ymin:ymax,xmin:xmax) = bwdist(mask(ymin:ymax,xmin:xmax)) <= margin;
        
        % The cc is assigned a new label, unless it has already been merged 
        if cc(r).labels(i) == 0
            cc(r).labels(i) = maxLabel;
        end
        
        % Find groups at smaller scales that can potentially be merged 
        mergedLabels = cc(r).labels(i);
        for rr=(r-1):-1:max(1,r-3)
            for j=1:cc(rr).NumObjects
                if ~any(mergedLabels == cc(rr).labels(j)) && any(mask(cc(rr).PixelIdxList{j}))
                    mergedLabels = [mergedLabels, cc(rr).labels(j)];
                end
            end
        end    
        
        % Merge labels (use the smallest label as the common label)
        commonLabel = min(mergedLabels);
        for rr=1:r
            cc(rr).labels(ismember(cc(rr).labels, mergedLabels)) = commonLabel;
        end
                        
        % Merge ccs at the same scale 
        for j=(i+1):cc(r).NumObjects
            if any(mask(cc(r).PixelIdxList{j}))
                cc(r).labels(j) = cc(r).labels(i);
            end
        end
        
        % If the component has not been merged, increase maxLabel
        if cc(r).labels(i) == maxLabel
            maxLabel = maxLabel + 1;
        end
    end
end

% Construct label map
branches = zeros(H,W,'uint8');
for r=1:R
    for i=1:cc(r).NumObjects
        branches(cc(r).PixelIdxList{i}) = cc(r).labels(i);
    end
end

% Adjust labels
oldLabels = unique(cat(2, cc(:).labels));
newLabels = 0:numel(oldLabels);
for i=2:numel(oldLabels) % 1st label is always zero (background)
    branches(branches == oldLabels(i)) = newLabels(i);
end