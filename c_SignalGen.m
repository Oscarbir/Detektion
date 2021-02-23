classdef c_SignalGen
    
    
    properties
        centerFq;
        bw;
        samples;
    end
    
    methods
        function obj = c_SignalGen(centerFq,bw)
            obj.centerFq=centerFq;
            obj.bw=bw;
            
        end
        
        function out = generateSignal(obj)
            
        end
    end
end

