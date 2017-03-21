classdef AMAT < handle
    properties
        scales = 2:41  
        ws = 1e-4      
        vistop = 0    
        filters
    end
    
    methods
        function obj = AMAT(arg)
%             if nargin > 0
%                 obj = copyFromObject(obj,arg);
%             else
%                 obj = initializeObject();
%             end
        end
    end
    
    methods(Access=private)
    end
end