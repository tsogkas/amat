classdef ScaleLoss < dagnn.Loss
% Most methods and properties are inherited from dagnn.Loss.
% We only need to create the scale-associated versions of the groundtruth
% labels and compute the respective instanceWeight arrays on the fly.
    properties
        scale = 0;
    end

    methods
        function outputs = forward(obj, inputs, params)
            % inputs{2} is the groundtruth HxWx1xB label map, where B is
            % the batchsize and inputs{2} values range from 0 (ignored
            % pixels for loss computation) to obj.numLabels with
            % obj.numLabels being the background label.
            % Create scale-associated versions of the groundtruth, that 
            % will be used to compute the scale-associated losses in the 
            % DSN layers by setting labels > s to the "background" label.
            [H,W,C,B] = size(inputs{1});
            inputs{2}(inputs{2}>obj.scale) = 1;
            obj.opts{end+1} = 'instanceWeights'; 
            obj.opts{end+1} = instanceWeights(inputs{2},obj.scale)/(H*W);
            outputs{1} = vl_nnloss(inputs{1}, inputs{2}, [], 'loss', obj.loss, obj.opts{:}) ;
            n = obj.numAveraged ;
            m = n + size(inputs{1},4) ;
            obj.average = (n * obj.average + gather(outputs{1})) / m ;
            obj.numAveraged = m ;
        end
        
        function [derInputs, derParams] = backward(obj, inputs, params, derOutputs)
            inputs{2}(inputs{2}>obj.scale) = 1;
            derInputs{1} = vl_nnloss(inputs{1}, inputs{2}, derOutputs{1}, 'loss', obj.loss, obj.opts{:}) ;
            derInputs{2} = [] ;
            derParams = {} ;
        end
        
        function obj = ScaleLoss(varargin)
            obj.load(varargin) ;
        end
    end
end


function w = instanceWeights(labels, numLabels)
% Returns a HxWx1xB matrix of pointwise weights according to the scheme
% suggested in the paper "Object Skeleton Extraction in Natural Images by
% Fusing Scale-associated Deep Side Outputs".
B = size(labels,4);
w = zeros(size(labels), 'like',labels); % HxWx1xB
b = zeros(numLabels,B,'like',labels);
for l=1:numLabels
    b(l,:) = sum(sum(labels==l));
end
assert(all(sum(b) == size(labels,1)*size(labels,2)))
b = 1./b; b(isinf(b)) = 0;
b = bsxfun(@rdivide, b, sum(b));
for l=1:numLabels
    w = w + bsxfun(@times,labels == l,reshape(b(l,:),1,1,1,[]));
end
end