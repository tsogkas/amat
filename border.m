function b = border(img, width)
% return a mask that corresponds to the border of the input. 
% img: 2D or 3D array (grayscale or color image)
% width: width of the border from all sides

if isscalar(img)
    error('Input must be an array or a vector of dimensions')
elseif isvector(img)
    sz = img;
else
    sz = size(img);
end
rows = sz(1);
cols = sz(2);
rest = prod(sz(3:end));
b = false(rows,cols,rest);
b(1:width,:,:) = true;
b(:,1:width,:) = true;
b((end-width+1):end,:,:) = true;
b(:,(end-width+1):end,:) = true;
b = reshape(b,sz);
    