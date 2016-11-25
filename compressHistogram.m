function  hout = compressHistogram(hin,K)
% COMPRESSHISTOGRAM Compresses a histogram encoding by keeping only the top
%   K bins with highest probability and setting all the rest to zero. 
%   NOTE: COMPRESSHISTOGRAM does not change the number of bins of the
%   histogram.
% 
%   hout = COMPRESSHISTOGRAM(hin,K) takes as input a HxWxCxBxR histogram 
% 	and returns its "compressed" version hout with the same dimensions, but
% 	with some of its bins zeroed out. If K is a real number in [0,1], then
% 	the top-K*100 % of the bins with highest probabilities are kept,
% 	whereas if K is an integer in [1,B], the top-K bins are kept. The
% 	probabilities for the remaining bins are set to zero.
% 
%   See also: histogram
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

hout = hin;
[H,W,C,B,R] = size(hin);
[~,sortedBinIndices] = sort(hin,4,'descend');
if K >= 0 && K <= 1 
    [~,top] = max(cumsum((1:B)/B) >= K);
    hout(sortedBinIndices(:,:,:,top+1:end,:)) = 0;
elseif K >= 1 && K <= B
    hout(sortedBinIndices(:,:,:,K+1:end,:)) = 0;
else
    error(['K should either be a real number in [0,1]'...
        'or an integer in [1,B], where B is the number of histogram bins'])
end
