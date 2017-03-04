function e = seg2edges(seg,thresh)
% Input is the output of mat2seg().

if nargin > 1 && isscalar(thresh) 
    e = edge(seg,'canny',thresh);
else  
    e = zeros(size(seg));
    numSegments = max(seg(:));
    % Add edges of "foreground" (explicitly selected) segments
    for i=1:numSegments
        e(bwperim(seg==i)) = i;
    end    
    % Add edges of "background" (remaining) segments
    bg = bwlabel(seg == 0);
    numNewSegments = max(bg(:));
    for i=1:numNewSegments
        e(bwperim(bg == i)) = numSegments + i;
    end
    e(border(e,2)) = 0;
end