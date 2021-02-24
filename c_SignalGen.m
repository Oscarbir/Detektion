classdef c_SignalGen
    
    
    properties
        centerFq;
        bw;
        gain;
        fs;
        N;
        
    end
    
    methods
        function obj = c_SignalGen(sParam)
            obj.centerFq=sParam.centerFq;
            obj.bw=sParam.bw;
            obj.gain=sParam.gain;
            obj.fs=sParam.fs;
            obj.N=obj.fs;
        end
        
        function [fAxis,out] = generateSignal(obj)
            
            df=obj.N/obj.fs; % delta f of sample
            k=obj.bw/df; % number of smaples for bw
            
            noiseFloor=randn(obj.N,1)+1i*rand(obj.N,1); %noisefloor
            noiseSignal=randn(k,1)+ 1i*randn(k,1);
            
            %zeropad to match length of noisefloor, gain added to signal
            zPaddedNoiseSignal=[zeros(obj.N/2-k/2,1);noiseSignal+obj.gain;zeros(obj.N/2-k/2,1,1)];
            
            out=zPaddedNoiseSignal+noiseFloor;
            
            fAxis = (-obj.N/2:obj.N/2-1)*(obj.fs/obj.N);
        end
    end
end

