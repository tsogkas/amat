% -------------------------------------------------------------------------
function O = edgeOrient( E, r )
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
