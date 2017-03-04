function e = depth2edges(mat,edgeDepth)
e = mat.depth;
e(border(e,3)) = inf; % avoid responses near the image border
e = e <= edgeDepth;

