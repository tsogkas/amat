function f = circleFilterBank(r)
% CIRCLEFILTERBANK Build a filterbank of circular filters.
% 
%   f = CIRCLEFILTERBANK(r) takes as input r, which can be either a
%   scalar or a vector returns a cell array f, containing circular masks of
%   respective radii.
% 
% See also: CIRCLE
% 
% Stavros Tsogkas, <tsogkas@cs.toronto.edu>
% Last update: October 2016 

numScales = numel(r);
f = cell(numScales,1);
for i=1:numScales
    f{i} = circle(r(i));
end
