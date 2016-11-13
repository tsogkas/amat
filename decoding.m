function g = decodeImage(f,method)
% DECODEIMAGE Compute decodings g at every pixel in the image. 
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if ndims(f) == 5        % HxWxCxBxS (compress histogram to scalar and expand)
elseif ndims(f) == 4    % HxWxCxS (expand scalar to disk area)
else error('The encodings must be 4-D or 5-D arrays')
end
