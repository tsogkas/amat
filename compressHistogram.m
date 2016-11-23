function  hout = compressHistogram(hin)
% Takes as input a histogram encoding of an image, hin, with dimensions
% HxWxCxBxR, and returns its "compressed" version with dimensions
% HxWxCxKxR, where K < B.