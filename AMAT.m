classdef AMAT < handle
    % TODO: add private flags for profiling
    properties
        scales  = 2:41  
        ws      = 1e-4      
        vistop  = 0  
        shape   = 'disk'
        filters 
        input
        reconstruction
        encoding
        axes
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
                % Only add these when they are efficient enough
                %mat.group();
                %mat.simplify();
            end
        end
        
        
        function compute(mat)
            mat.computeEncodings();
            mat.computeCosts();
            mat.setCover();
        end
        
        function group(mat)
        end
        
        function simplify(mat)
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
            mat.axes           = zeroLabNormalized;
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
            mat.axes  = labNormalized2rgb(mat.axes);
            mat.reconstruction = mat2reconstruction(mat.axes,mat.radius,mat.depth,mat.scales);
            
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
                subplot(223); imshow(mat.axes); title('A-MAT axes')
                subplot(224); imshow(mat.radius,[]); title('A-MAT radii')
                drawnow;
            end
            
        end
        
        function visualize(mat)
            subplot(221); imshow(mat.axes);             title('Medial axes');
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
            % Fast version of imageEncoding, using convolutions with 
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