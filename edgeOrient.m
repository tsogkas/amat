function O = edgeOrient( E, r, matlabOnly)
% Actually the code for both subfunctions is taken from Piotr Dollar's edge
% detection code. The second one is the variant with matlab-only dependency

if nargin < 2, r = 4; end
if nargin < 3, matlabOnly = false; end
if matlabOnly
    O = edgeOrientMatlab(E,r);
else
    O = edgeOrientDollar(E,r);
end

% -------------------------------------------------------------------------
function O = edgeOrientDollar(E,r)
% -------------------------------------------------------------------------
[Ox,Oy]=gradient2(convTri(E,r));
[Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);

% -------------------------------------------------------------------------
function O = edgeOrientMatlab(E,r)
% -------------------------------------------------------------------------
% Compute very approximate orientation map from edge map. We have replaced
% convTri with the slightly slower equivalent code to reduce dependencies.
% E2=convTri(E,r); 
% ----------------------------------------------
if(r<=1), p=12/r/(r+2)-2; f=[1 p 1]/(2+p); r=1;
else f=[1:r r+1 r:-1:1]/(r+1)^2; end
E2 = padarray(E,[r r],'symmetric','both');
E2 = convn(convn(E2,f,'valid'),f','valid');
% ----------------------------------------------
f=[-1 2 -1];
Dx=conv2(E2,f,'same'); Dy=conv2(E2,f','same');
F=conv2(E2,[1 0 -1; 0 0 0; -1 0 1],'same')>0;
Dy(F)=-Dy(F); O=mod(atan2(Dy,Dx),pi);
