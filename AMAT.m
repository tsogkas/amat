classdef AMAT < handle
    % TODO: add private flags for profiling
    % TODO: double-check all result types (keep 'single' or switch to double?)
    properties
        scales  = 2:41  
        ws      = 1e-4      
        vistop  = 0  
        shape   = 'disk'
        filters 
        input
        reconstruction
        encoding
        axis
        radius
        depth
        price
        cost
        branches
        simplified
        covered
    end
    
    properties(Access=private)
    end
    
    methods
        function mat = AMAT(img,varargin)
            if nargin > 0
                assert(size(img,3)>=2, 'Input image must be 2D or 3D array')
                mat.initialize(img,varargin{:});
                mat.compute();
            end
        end
        
        
        function compute(mat)
            mat.computeEncodings();
            mat.computeCosts();
            mat.setCover();
        end
        
        function group(mat,marginFactor,colortol)
            if nargin < 2, marginFactor = 1; end
            if nargin < 3, colortol = 0.05; end
            
            % Compute individual radius maps and connected components
            R = numel(mat.scales);
            for r=R:-1:1
                cc(r) = bwconncomp(mat.radius == mat.scales(r));
            end
            
            % Initialize mask and maxLabel
            [H,W,C] = size(mat.input);
            mask = false(H,W);  % proximity mask
            maxLabel = 1;       % initialize maxLabel
            % Convert to Lab and reshape axis encodings for convenience
            if colortol
                mataxis = reshape(rgb2labNormalized(mat.axis), H*W,C);
            end
            
            % For all scales
            for r=1:R
                cc(r).labels = zeros(1, cc(r).NumObjects); % zero for non-examined ccs
                margin = ceil(marginFactor*r)+1;
                % For all connected components at the same scale
                for i=1:cc(r).NumObjects;
                    % Create proximity mask in rectangle around cc for efficiency
                    mask(:) = false; mask(cc(r).PixelIdxList{i}) = true;
                    idxcc = cc(r).PixelIdxList{i};
                    [y,x] = ind2sub([H,W], idxcc);
                    xmin = max(1,min(x)-margin); xmax = min(W,max(x)+margin);
                    ymin = max(1,min(y)-margin); ymax = min(H,max(y)+margin);
                    mask(ymin:ymax,xmin:xmax) = bwdist(mask(ymin:ymax,xmin:xmax)) <= margin;
                    
                    % The cc is assigned a new label, unless it has already been merged
                    if cc(r).labels(i) == 0
                        cc(r).labels(i) = maxLabel;
                    end
                    
                    % Find groups at smaller scales that can potentially be merged
                    mergedLabels = cc(r).labels(i);
                    for rr=(r-1):-1:max(1,r-4)
                        for j=1:cc(rr).NumObjects
                            if ~any(mergedLabels == cc(rr).labels(j)) && merge(cc(rr).PixelIdxList{j})
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
                        if merge(cc(r).PixelIdxList{j})
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
            matbranches = zeros(H,W);
            for r=1:R
                for i=1:cc(r).NumObjects
                    matbranches(cc(r).PixelIdxList{i}) = cc(r).labels(i);
                end
            end
            
            % Adjust labels. We do not need to explicitly remove the zero labels
            % because cc.labels() does not include any zero (0) labels.
            oldLabels = unique(cat(2, cc(:).labels));
            newLabels = 1:numel(oldLabels);
            for i=1:numel(oldLabels)
                matbranches(matbranches == oldLabels(i)) = newLabels(i);
            end
            mat.branches = matbranches;
            
            % Nested functions --------------------------------------------
            function res = merge(idx)
                res = isCloseSpace(idx);
                if colortol, res = res && isCloseColor(idx); end
            end
            
            function res = isCloseSpace(idx)
                res = any(mask(idx));
            end
            
            function res = isCloseColor(idx)
                res = norm( mean(mataxis(idxcc,:),1)-...
                    mean(mataxis(idx,:),1) ) < colortol;
            end
            
            
        end
        
        function simplify(mat,method,param)
            % Default input arguments
            if nargin < 3, param  = 3; end
            if nargin < 2, method = 'dilation'; end
            
            % Post-processing function
            switch method
                case 'dilation'
                    SE = strel('disk',param);
                    process = @(x) imdilate(x,SE);
                case 'iso-dilation'
                    process = @(x) bwdist(x) <= param;
                case 'skeleton'
                    process = @(x) bwmorph(x,'skel',inf);
                case 'afmm-skeleton'
                    process = @(x) skeleton(x)>=param;
                otherwise
                    error(['Method not supported. Supported methods are:\n' ...
                        'dilation, iso-dilation, skeleton, afmm-skeleton.'])
            end
            
            % The group labels are already sorted and first label is zero (background)
            numBranches = max(mat.branches(:));
            [H,W,C]     = size(mat.input);
            matbranches = zeros(H,W);
            matradius   = zeros(H,W);
            for i=1:numBranches
                % Old branch points, radii, and respective cover.
                branchOld = mat.branches == i;
                radiusOld = branchOld .* double(mat.radius);
                cover     = mat2mask(radiusOld, mat.scales)>0;
                % Apply post-processing and thinning to selected branch.
                % Crop the necessary region for more efficiency.
                % Dilation and iso-dilation are applied on the branch points, whereas
                % skeletonization is applied on the cover mask.
                if strcmp(method,'dilation') || strcmp(method,'iso-dilation')
                    branchNew = bwmorph(process(branchOld),'thin',inf);
                else
                    branchNew = bwmorph(process(cover),'thin',inf);
                end
                
                % Compute new radii as distance transform on reconstructed cover.
                radiusNew = bwdist(bwperim(cover)).* double(branchNew);
                % Find closest radii in the subset of the acceptable scale values.
                valid = radiusNew > 0;
                [~,idx] = min(abs(bsxfun(@minus,radiusNew(valid),mat.scales)),[],2);
                radiusNew(valid) = mat.scales(idx);
                % Assign values in global label and radius map.
                matbranches(valid) = i;
                matradius(valid) = radiusNew(valid);
            end
            assert(all(matbranches(matbranches>0) & matradius(matbranches>0)))
            assert(all(ismember(matradius(matradius>0), mat.scales)))
            
            % Make sure there are no gaps among branch labels
            newLabels = unique(matbranches); newLabels(1) = []; % first group is zero
            for i=1:numel(newLabels)
                matbranches(matbranches == newLabels(i)) = i;
            end
            
            % Find which pixels have been removed and which have been added
            oldpts  = any(mat.axis,3);
            newpts  = matbranches > 0;
            removed = oldpts & ~newpts;
            added   = newpts & ~oldpts;
            
            % Update depth
            % NOTE: there is a discrepancy between
            % mat2mask(double(newpts).*radius,mat.scales) and
            % mat.depth + depthAdded - depthRemoved. This is probably because when the
            % new radii are changed EVEN FOR THE POINTS THAT ARE NOT REMOVED.
            % depthAdded   = mat2mask(radius.*double(added),       mat.scales);
            % depthRemoved = mat2mask(mat.radius.*double(removed), mat.scales);
            matdepth = mat2mask(matradius,mat.scales);
            
            % Update MAT encodings
            [y,x] = find(newpts);
            r   = matradius(newpts);
            R   = numel(mat.scales);
            enc = reshape(permute(imageEncoding(rgb2labNormalized(mat.input),mat.scales),[1 2 4 3]), [], C);
            rind= containers.Map(mat.scales,1:numel(mat.scales));
            for i=1:numel(r), r(i) = rind(r(i)); end % map scales to scale indexes
            idx = sub2ind([H,W,R], y(:),x(:),r(:));
            newaxis = reshape(rgb2labNormalized(zeros(H,W,C)),H*W,C);
            newaxis(newpts,:) = enc(idx,:); % remember that encodings are in LAB!
            newaxis = labNormalized2rgb(reshape(newaxis,H,W,C));
            
            % Update reconstruction
            matreconstruction = mat2reconstruction(reshape(newaxis,H,W,C),...
                matradius, matdepth, mat.scales);
            
            % Update mat fields
            mat.radius   = matradius;
            mat.branches = matbranches;
            mat.axis     = newaxis;
            mat.depth    = matdepth;
            mat.reconstruction = matreconstruction;

        end
        
        function computeEncodings(mat)
            switch mat.shape
                case 'disk'
                    computeDiskEncodings(mat);
                otherwise, error('Invalid shape')
            end            
        end
        
        function computeCosts(mat)
            switch mat.shape
                case 'disk'
                    computeDiskCosts(mat);
                otherwise, error('Invalid shape')
            end
        end
        
        function setCover(mat)
            % -------------------------------------------------------------
            % Greedy approximation of the weighted set cover problem.
            % -------------------------------------------------------------
            % - Disk cost: cost incurred by selecting ena r-disk, centered at (i,j).
            % - numNewPixelsCovered: number of NEW pixels covered by a selected disk.
            % - Cost per pixel: diskCost / numNewPixelsCovered.
            % - Disk cost effective: adjusted normalized cost: costPerPixel + scaleTerm
            %       where scaleTerm is a term that favors selecting disks of larger
            %       radii. Such a term is necessary, to resolve selection of disks in
            %       the case where diskCost is zero for more than on radii.
            %
            % TODO: is there a way to first sort scores and then pick the next one in
            %       queue, to avoid min(diskCostEffective(:)) in each iteration?
            
            % Initializations
            [H,W,C,R]          = size(mat.encoding);
            zeroLabNormalized  = rgb2labNormalized(zeros(H,W,C,'single'));
            mat.input          = reshape(img, H*W, C);
            mat.reconstruction = reshape(zeroLabNormalized,H*W,C);
            mat.axis           = zeroLabNormalized;
            mat.radius         = zeros(H,W,'single');
            mat.depth          = zeros(H,W,'single'); % #disks points(x,y) is covered by
            mat.price          = zeros(H,W,'single'); % error contributed by each point
            mat.covered        = false(H,W);
            % Flag border pixels that cannot be accessed by filters.
            r = mat.scales(1);
            mat.covered([1:r,end-r+1:end], [1,end]) = true;
            mat.covered([1,end], [1:r,end-r+1:end]) = true;
            BIG = 1e30;
            
            % Compute how many pixels are covered be each r-disk.
            diskAreas = cellfun(@nnz,mat.filters);
            diskCost  = mat.cost;
            numNewPixelsCovered = repmat(reshape(diskAreas,1,1,[]), [H,W]);
            
            % Add scale-dependent cost term to favor the selection of larger disks.
            costPerPixel = diskCost ./ numNewPixelsCovered;
            diskCostEffective = bsxfun(@plus, costPerPixel, ...
                reshape(mat.ws ./ mat.scales, 1,1,[]));
            
            % Print remaining pixels to be covered in these points
            printBreakPoints = floor((4:-1:1).*(H*W/5));
            
            fprintf('Pixels remaining: ');
            [x,y] = meshgrid(1:W,1:H);
            while ~all(mat.covered(:))
                % Find the most cost-effective disk at the current iteration
                [minCost, indMin] = min(diskCostEffective(:));
                if isinf(minCost),
                    warning('Stopping: selected disk has infinite cost.')
                    break;
                end
                
                [yc,xc,rc] = ind2sub([H,W,R], indMin);
                D = (x-xc).^2 + (y-yc).^2 <= mat.scales(rc)^2; % points covered by the selected disk
                newPixelsCovered  = D & ~mat.covered;      % NEW pixels that are covered by D
                if ~any(newPixelsCovered(:))
                    warning('Stopping: selected disk covers zero (0) new pixels.')
                    break;
                end
                
                % Update MAT
                update(mat)
                % Update costs
                [yy,xx] = find(newPixelsCovered);
                xmin = min(xx); xmax = max(xx);
                ymin = min(yy); ymax = max(yy);
                newPixelsCovered = double(newPixelsCovered);
                for r=1:R
                    scale = mat.scales(r);
                    x1 = max(xmin-scale,1); y1 = max(ymin-scale,1);
                    x2 = min(xmax+scale,W); y2 = min(ymax+scale,H);
                    % Find how many of the newPixelsCovered are covered by other disks.
                    numPixelsSubtracted = ...
                        conv2(newPixelsCovered(y1:y2,x1:x2), mat.filters{r},'same');
                    % and subtract the respective counts from those disks.
                    numNewPixelsCovered(y1:y2,x1:x2, r) = ...
                        numNewPixelsCovered(y1:y2,x1:x2, r) - numPixelsSubtracted;
                    % update diskCost, costPerPixel, and diskCostEfficiency *only* for
                    % the locations that have been affected, for efficiency.
                    diskCost(y1:y2,x1:x2, r) = diskCost(y1:y2,x1:x2, r) - ...
                        numPixelsSubtracted .* costPerPixel(y1:y2,x1:x2, r);
                    costPerPixel(y1:y2,x1:x2, r) = diskCost(y1:y2,x1:x2, r) ./ ...
                        max(eps,numNewPixelsCovered(y1:y2,x1:x2, r)) + ... % avoid 0/0
                        BIG*(numNewPixelsCovered(y1:y2,x1:x2, r) == 0);    % x/0 = inf
                    diskCostEffective(y1:y2,x1:x2, r) = ...
                        costPerPixel(y1:y2,x1:x2, r) + mat.ws/mat.scales(r);
                end
                % Make sure the same point is not selected again
                diskCost(yc,xc,:) = BIG; diskCostEffective(yc,xc,:) = BIG;
                
                
                if mat.vistop, visualizeProgress(mat,diskCostEffective); end
                if ~isempty(printBreakPoints) && nnz(~mat.covered) < printBreakPoints(1)
                    fprintf('%d...',printBreakPoints(1))
                    printBreakPoints(1) = [];
                end
            end
            fprintf('\n')
            mat.input = labNormalized2rgb(reshape(mat.input,H,W,C));
            mat.axis  = labNormalized2rgb(mat.axis);
            mat.reconstruction = mat2reconstruction(mat.axis,mat.radius,mat.depth,mat.scales);
            
            function update(mat)
                mat.price(newPixelsCovered) = minCost / numNewPixelsCovered(yc,xc,rc);
                mat.covered(newPixelsCovered) = true;
                mat.depth(D) = mat.depth(D) + 1;
                mat.axis(yc,xc,:) = mat.encoding(yc,xc,:,rc);
                mat.radius(yc,xc) = mat.scales(rc);    
            end
            function visualizeProgress(mat,diskCost)
                % Sort costs in ascending order to visualize updated top disks.
                [~, indSorted] = sort(diskCost(:),'ascend');
                [yy,xx,rr] = ind2sub([H,W,R], indSorted(1:mat.vistop));
                subplot(221); imshow(reshape(mat.input, H,W,[]));
                viscircles([xc,yc],rc, 'Color','k','EnhanceVisibility',false); title('Selected disk');
                subplot(222); imshow(bsxfun(@times, reshape(mat.input,H,W,[]), double(~mat.covered)));
                viscircles([xx,yy],rr,'Color','w','EnhanceVisibility',false,'Linewidth',0.5);
                viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false);
                viscircles([xc,yc],rc,'Color','y','EnhanceVisibility',false);
                title(sprintf('K: covered %d/%d, W: Top-%d disks,\nB: Top-1 disk, Y: previous disk',...
                    nnz(mat.covered),H*W,mat.vistop))
                subplot(223); imshow(mat.axis); title('A-MAT axes')
                subplot(224); imshow(mat.radius,[]); title('A-MAT radii')
                drawnow;
            end
            
        end
        
        function visualize(mat)
            subplot(221); imshow(mat.axis);             title('Medial axes');
            subplot(222); imshow(mat.radius,[]);        title('Radii');
            subplot(223); imshow(mat.input);            title('Original image');
            subplot(224); imshow(mat.reconstruction);   title('Reconstructed image');
        end
    end
    
    methods(Access=private)
        function initialize(mat,img,opts)
            defaults = {'scales',   2:41,...
                        'ws',       1e-4,...
                        'vistop',   0,...
                        'shape',    'disk'
                        };
            opts = parseVarargin(defaults,opts);
            if isscalar(opts('scales'))
                mat.scales  = 2:opts('scales');
            else
                mat.scales  = opts('scales');
            end
            mat.ws      = opts('ws');
            mat.vistop  = opts('vistop');
            mat.shape   = opts('shape');
            mat.input   = rgb2labNormalized(im2double(img));
            mat.initializeFilters();
            
            
        end
        
        function initializeFilters(mat)
            numScales = numel(mat.scales);
            mat.filters = cell(1, numScales);
            switch mat.shape
                case 'disk'
                    f = @(x) disk(x);
                otherwise, error('Invalid filter shape')
            end
            for i=1:numScales
                mat.filters{i} = double(f(mat.scales(i))); 
            end
        end
        
        function computeDiskEncodings(mat)
            % Efficient implementation, using convolutions with 
            % circles + cumsum instead of convolutions with disks.
            [H,W,C] = size(mat.input); R = numel(mat.scales);
            cfilt = cell(1,R); cfilt{1} = double(disk(mat.scales(1)));
            for r=2:R, cfilt{r} = double(circle(mat.scales(r))); end
            enc = zeros(H,W,C,R);
            for c=1:C
                for r=1:R
                    enc(:,:,c,r) = conv2(mat.input(:,:,c),cfilt{r},'same');
                end
            end
            enc   = cumsum(enc,4);
            areas = cumsum(cellfun(@nnz,cfilt));
            enc   = bsxfun(@rdivide, enc, reshape(areas,1,1,1,[]));
            mat.encoding = enc;
        end
        
        function computeDiskCosts(mat)
            % This function computes a heuristic that represents the 
            % ability to reconstruct a disk-shaped part of the input image
            % using the mean RGB values computed over the same area.
            % Intuitively, the idea behind this heuristic is the following:
            % In order to accurately reconstruct an image disk of radius r 
            % using its mean RGB values, we must also be able to reconstruct 
            % *every* fully contained disk of radius r' < r 
            % (uniformity criterion).
            %
            % Terminology: an r-disk is a disk of radius = r.
            %
            % The heuristic we use sums all square errors between the 
            % encoding of an r-disk centered at a point (i,j) and the 
            % encodings of all FULLY CONTAINED disks. 
            % Written in a simplified mathematical form, for a given r_k-disk:
            % M_rk = sum_i(I_rk)/D_rk; M_ri = sum_i(I_ri)/R_ri;
            % Cost = sum((M_rk-M_ri)^2) for all enclosed ri-disks.
            % Cost = sum( M_rk^2 + M_ri^2 - 2*M_rk*M_ri ) = ...
            % D_rk*enc2 + conv2(enc2) + 2 .* enc .* conv2(enc)
            % Given an r-disk, filters(r-i+1) is a mask that marks the 
            % centers of all contained i-disks.
            [H,W,C,R] = size(mat.encoding);
            cfilt     = cell(1,R);
            cfilt{1}  = double(disk(mat.scales(1)-1));
            for r=2:R, cfilt{r} = double(circle(mat.scales(r-1))); end
            % Precompute necessary quantitities. We use circular filters applied on
            % cumulative sums instead of disk filters, for efficiency.
            enc      = mat.encoding;
            enc2     = enc.^2;
            enccsum  = cumsum(enc,4);
            enc2csum = cumsum(enc2,4);
            nnzcd    = cumsum(cumsum(cellfun(@nnz, cfilt)));
            
            diskCost = zeros(H,W,C,R);
            for c=1:C
                for r=1:R
                    sumMri  = zeros(H,W);
                    sumMri2 = zeros(H,W);
                    for i=1:r
                        sumMri  = sumMri  + conv2(enccsum(:,:,c,i), cfilt{r-i+1},'same');
                        sumMri2 = sumMri2 + conv2(enc2csum(:,:,c,i),cfilt{r-i+1},'same');
                    end
                    diskCost(:,:,c,r) = enc2(:,:,c,r)*nnzcd(r) + sumMri2 - 2*enc(:,:,c,r).*sumMri;
                end
            end
            
            % Fix boundary conditions. Setting scale(r)-borders to a very big cost
            % helps us avoid selecting disks that cross the image boundaries.
            % We do not use Inf to avoid complications in the greedy set cover
            % algorithm, caused by inf-inf subtractions and inf/inf divisions.
            % Also, keep in mind that max(0,NaN) = 0.
            BIG = 1e30;
            for r=1:R
                scale = mat.scales(r);
                diskCost([1:scale, end-scale+1:end],:,:,r) = BIG;
                diskCost(:,[1:scale, end-scale+1:end],:,r) = BIG;
            end
            
            % Sometimes due to numerical errors, cost are slightly negative. Fix this.
            diskCost = max(0,diskCost);
            
            % Combine costs from different channels
            if C > 1
                wc = [0.5,0.25,0.25]; % weights for luminance and color channels
                diskCost = diskCost(:,:,1,:)*wc(1) + diskCost(:,:,2,:)*wc(2) + diskCost(:,:,3,:)*wc(3);
                diskCost = squeeze(diskCost);
            end
            mat.cost = diskCost;
        end
        
    end
    
end