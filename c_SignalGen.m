classdef c_SignalGen
    
    
    properties
        centerFq;
        bw;
        power;
        fs;
        N;
        noise;
        
    end
    
    methods
        function obj = c_SignalGen(sParam)
            obj.centerFq=sParam.centerFq;
            obj.bw=sParam.bw;
            obj.power=sParam.power;
            obj.fs=sParam.fs;
            obj.N=sParam.N;
            obj.noise=sParam.noise;
        end
        
        function [fAxis,out] = generateSignal(obj)
            
            
            df=obj.fs/obj.N; % delta f of sample
            k=obj.bw/df; % number of smaples for bw
            
            noise=randn(obj.N,1)+1i*rand(obj.N,1); %noisefloor
            randSignal=(randn(k,1)+ 1i*randn(k,1)); % random signal with bw, in fqdomain
            
            %padd signal and noise with zeros
            randSignalZP=[zeros(obj.N/2-k/2,1); randSignal  ; zeros(obj.N/2-k/2,1,1)];
            noiseZP=[noise(1:obj.N/2-k/2);zeros(k,1);noise(obj.N/2+k/2:obj.N-1)];
            
            %convert signals to timeDomain
            outSignal=(ifft(fftshift(randSignalZP)));
            outNoise=obj.noise*(ifft(fftshift(noiseZP)));
            
            % gain calc for power specified
            dt=1/obj.fs;
            
            %calc and normalize energy in signal to 1 w, for easy scaling
            E=sum((abs(outSignal).^2.*dt));
            outSignal_scaled=outSignal./sqrt(E);
            
            %scale signal for target power
            E=obj.power*length(outSignal_scaled)*dt;
            outSignal_scaled=outSignal_scaled.*sqrt(E);
            out=outSignal_scaled+outNoise;
            
            fAxis = (-obj.N/2:obj.N/2-1)*(obj.fs/obj.N);
        end
    end
end

