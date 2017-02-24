% %% Compute radius maps and connected components
% img = imresize(imread('google.jpg'),0.5,'bilinear');
% mat = amat(img);
% R = 40;
% [H,W,C] = size(mat.input);
% 
% rad = zeros([H,W,R]);
% cc = cell(1,R);
% for r=1:R
%     rad(:,:,r) = mat.radius == r;
%     cc{r} = bwconncomp(rad(:,:,r));
% end

%% Grouping scheme
branches = zeros(H,W,'uint8');
maxLabel = 1;
for r=1:10
    % Initialize grouping be assigning each conn component its own label
    cc{r}.labels = zeros(1, cc{r}.NumObjects);
    for i=1:cc{r}.NumObjects;
        % Indices of currently selected connected component
        [y,x] = ind2sub([H,W], cc{r}.PixelIdxList{i});
        % Create box area of specified margin around component
        margin = ceil(0.7*r);
        xmin = max(1,min(x)-margin); xmax = min(W,max(x)+margin);
        ymin = max(1,min(y)-margin); ymax = min(H,max(y)+margin);
        % By default the component is assigned a new label, unless it has
        % already been merged with another component
        if cc{r}.labels(i) == 0
            cc{r}.labels(i) = maxLabel;
            % First check if there are components of smaller radii nearby
            if r > 1
                for j=1:cc{r-1}.NumObjects
                    [yy,xx] = ind2sub([H,W], cc{r-1}.PixelIdxList{j});
                    xcontained = xx >= xmin & xx <= xmax;
                    ycontained = yy >= ymin & yy <= ymax;
                    % if the connected component is close enough, merge them.
                    % if there is a match with > 1 components @r-1, then these
                    % should have the same label anyway.
                    if any(xcontained & ycontained)
%                         tmp = zeros(H,W,'uint8'); 
%                         tmp(cc{r}.PixelIdxList{i}) = 1;
%                         tmp(cc{r-1}.PixelIdxList{j}) = 2;
%                         imagesc(tmp); axis off image; 
%                         title('To be merged @previous scale'); 
%                         rectangle('Position',[xmin,ymin,xmax-xmin+1,ymax-ymin+1],'EdgeColor','y')
%                         keyboard;                        
                        cc{r}.labels(i) = cc{r-1}.labels(j);
                    end
                end
            end
            % If the component has not been merged, increase maxLabel
            if cc{r}.labels(i) == maxLabel, maxLabel = maxLabel + 1; end
        end
        % Now check if components at the same scale are nearby
        for j=i+1:cc{r}.NumObjects;
            [yy,xx] = ind2sub([H,W], cc{r}.PixelIdxList{j});
            xcontained = xx >= xmin & xx <= xmax;
            ycontained = yy >= ymin & yy <= ymax;
            if any(xcontained & ycontained)
%                 tmp = zeros(H,W,'uint8');
%                 tmp(cc{r}.PixelIdxList{i}) = 1;
%                 tmp(cc{r}.PixelIdxList{j}) = 2;
%                 imagesc(tmp); axis off image; 
%                 title('To be merged @same scale');
%                 rectangle('Position',[xmin,ymin,xmax-xmin+1,ymax-ymin+1],'EdgeColor','y')
%                 keyboard; 
                cc{r}.labels(j) = cc{r}.labels(i);
            end
        end
    end
    % Update label map
    for i=1:cc{r}.NumObjects
        branches(cc{r}.PixelIdxList{i}) = cc{r}.labels(i);
    end
end
imagesc(branches)