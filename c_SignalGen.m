classdef c_SignalGen
    
    
    properties
        centerFq;
        bw;
        gain;
        fs;
        N;
        noiseFloor;
        
    end
    
    methods
        function obj = c_SignalGen(sParam)
            obj.centerFq=sParam.centerFq;
            obj.bw=sParam.bw;
            obj.gain=sParam.gain;
            obj.fs=sParam.fs;
            obj.N=sParam.N;
            obj.noiseFloor=sParam.noiseFloor;
        end
        
        function [fAxis,out] = generateSignal(obj)
            
            df=obj.fs/obj.N; % delta f of sample
            k=obj.bw/df; % number of smaples for bw
            
            noise=randn(obj.N,1)+1i*rand(obj.N,1); %noisefloor
            randSignal=randn(k,1)+ 1i*randn(k,1); % random signal with bw, in fqdomain
            
            %padd signal and noise with zeros
            randSignalZP=[zeros(obj.N/2-k/2,1); randSignal ; zeros(obj.N/2-k/2,1,1)];
            noiseZP=[noise(1:obj.N/2-k/2);zeros(k,1);noise(obj.N/2+k/2:obj.N-1)];
            
            %convert signals to timeDomain
            outSignal=obj.gain*(ifft(fftshift(randSignalZP)));
            outNoise=obj.noiseFloor*(ifft(fftshift(noiseZP)));
            
            out=outSignal+outNoise;
            
            fAxis = (-obj.N/2:obj.N/2-1)*(obj.fs/obj.N);
        end
    end
end

