function [rec, comp] = spbmil2reconstruction(img)
histf = computeHistogramFeatures(img);
spb = spbMIL(img);
% Create rectangular filters at all scales and orientations
filters = cell(numel(histf.thetas), numel(histf.scales));
for s=1:numel(histf.scales)
    sc = histf.scales(s);
    rect = ones(2*sc+1, 2*histf.opts.ratio*sc+1);
    for o=1:numel(histf.thetas)
        filters{o,s} = imrotate(rect,rad2deg(histf.thetas(o)));
    end
end
% Pad input image, scales and thetas by replicating image border. This is
% necessary to compute mean value encoding correctly, and to avoid out of
% bound errors. Pad = 2*max(scales) to account for rotated versions of filters.
img    = im2double(img);
pad    = round(max(vec(cellfun(@(x)max(size(x)), filters)))/2);
imgPad = padarray(img,[pad,pad],'replicate','both');
pb     = padarray(spb.thin,[pad,pad],0,'both');
scales = padarray(spb.scalesMap,[pad,pad],1,'both');
thetas = padarray(spb.orientMap,[pad,pad],1,'both');
% Sort medial point locations, scales and orientations by their scores.
[pbSorted, indSorted] = sort(pb(:),'descend');
scales = scales(indSorted);
thetas = thetas(indSorted);

% Selecting medial points with the highest scores, compute mean values over
% the respective filter support and stop when the entire image has been covered.
[H,W,C] = size(imgPad);
rec = zeros(H,W,C);
covered = zeros(H,W); covered(border(covered,pad)) = 1;
[y,x] = ind2sub([H,W], indSorted);
i = 1;
while pbSorted(i) > 0.1 && ~all(covered(:))
% while pbSorted(i) > 0.1
    yy = y(i); xx = x(i);
    sc = scales(i);
    th = thetas(i);
    rf = filters{th,sc};
    [h,w] = size(rf);
    hs = floor(h/2); y1 = yy-hs; y2 = y1+h-1;
    ws = floor(w/2); x1 = xx-ws; x2 = x1+w-1;
    assert(y2-y1+1 == h);
    assert(x2-x1+1 == w);
    % Compute mean value encoding of local rectangular patch
    patch = bsxfun(@times, imgPad(y1:y2,x1:x2,:), rf);
    mval = sum(sum(patch))/nnz(rf);
    % Only fill in pixels that have no<t been already covered
    patch = bsxfun(@times, mval, rf );
    rec(y1:y2,x1:x2,:) = rec(y1:y2,x1:x2,:)+patch;
    covered(y1:y2,x1:x2) = covered(y1:y2,x1:x2) + rf;
    i = i+1;
end
rec = bsxfun(@rdivide, rec, covered);
rec = rec(pad+1:end-pad,pad+1:end-pad,:);
% Fill in holes that might be left
rec = reshape(inpaint_nans(rec), size(img,1), size(img,2),3);
comp = i / (size(img,1)*size(img,2));
