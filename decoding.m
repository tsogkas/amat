function g = decoding(f,method)
% DECODING Compute decodings g across the image for all filters. 
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

if ndims(f) == 5
elseif ndims(f) == 4
else error('The encodings must be 4-D or 5-D arrays')
end
