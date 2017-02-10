classdef Slice < dagnn.ElementWise
    properties
        dim = 3
        splitPoints = []
    end
    
    properties (Transient)
        inputSizes  = {};
    end
    
    methods
        function outputs = forward(obj, inputs, params)
            % by default split each block on its own
            if isempty(obj.splitPoints)
                obj.splitPoints = 1:size(inputs{1},obj.dim);
            end
            outputs = vl_nnslice(inputs, obj.dim, obj.splitPoints) ;
            obj.inputSizes = cellfun(@size, inputs, 'UniformOutput', false) ;
        end
        
        function [derInputs, derParams] = backward(obj, inputs, params, derOutputs)
            derInputs = vl_nnslice(inputs, obj.dim, obj.splitPoints, derOutputs{1}) ;
            derParams = {} ;
        end
        
        function reset(obj)
            obj.splitPoints = {} ;
            obj.inputSizes  = {} ;
        end
        
        function outputSizes = getOutputSizes(obj, inputSizes)
            sz = inputSizes{1};
            outputSizes = cell(1, numel(obj.splitPoints)+1);
            sz(obj.dim) = obj.splitPoints(1)-1; 
            outputSizes{1} = sz;
            for k = 2:numel(outputSizes)
                sz(obj.dim) = obj.splitPoints(k)-obj.splitPoints(k-1)+1;
                outputSizes{k} = sz;
            end
            outputSizes{end} = size(inputSizes,obj.dim)-obj.splitPoints(k)+1;
        end
        
        function rfs = getReceptiveFields(obj)
            numInputs = numel(obj.net.layers(obj.layerIndex).inputs) ;
            if obj.dim == 3 || obj.dim == 4
                rfs = getReceptiveFields@dagnn.ElementWise(obj) ;
                rfs = repmat(rfs, numInputs, 1) ;
            else
                for i = 1:numInputs
                    rfs(i,1).size = [NaN NaN] ;
                    rfs(i,1).stride = [NaN NaN] ;
                    rfs(i,1).offset = [NaN NaN] ;
                end
            end
        end
        
        function load(obj, varargin)
            s = dagnn.Layer.argsToStruct(varargin{:}) ;
            % backward file compatibility
            if isfield(s, 'numInputs'), s = rmfield(s, 'numInputs') ; end
            load@dagnn.Layer(obj, s) ;
        end
        
        function obj = Slice(varargin)
            obj.load(varargin{:}) ;
        end
    end
end


function y = vl_nnslice(inputs, dim, splitPoints, dzdy)
%VL_NNSLICE CNN split input into multiple outputs.
%  Y = VL_NNSPLIT(INPUTS, DIM, SPLITPOINTS) splits the inputs in the cell
%  array INPUTS along dimension DIM at split points SPLITPOINTS, generating
%  numel(SPLITPOINTS)+1 outputs.
%
%  DZDINPUTS = VL_NNSPLIT(INPUTS, DIM, SPLITPOINTS, DZDY) concatenates the
%  derivatives of the split blocks, contained in DZDY.
%
% Copyright (C) 2017 Stavros Tsogkas
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

if nargin < 2, dim = 3; end;
if nargin < 3, splitPoints = 1; end
if nargin < 4, dzdy = []; end;

if isempty(dzdy) % forward pass
    assert(numel(splitPoints) < size(inputs{1},dim))
    y = cell(1, numel(splitPoints)+1);
    s.type = '()' ;
    s.subs = {':', ':', ':', ':'} ;
    start = 1;
    for i=1:numel(splitPoints)
        stop = start+splitPoints(i);
        s.subs{dim} = start:stop-1;
        y{i} = subsref(inputs{1},s);
        start = stop;
    end
    % Last blob
    s.subs{dim} = start:size(inputs{1},dim);
    y{end} = subsref(inputs{1},s);
else
    y = cat(dim, dzdy{:});
end
end
