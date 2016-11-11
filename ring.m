function r = ring(rin, rout)
[x,y] = meshgrid(-rout:rout, -rout:rout);
r = (x.^2 + y.^2 > rin^2) & (x.^2 + y.^2 <= rout^2) ;