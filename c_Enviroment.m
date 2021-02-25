classdef c_Enviroment
    %coordsys(0,0) at reciever, distance in m. (x,y). 
    % todo implement noisefloor
    properties
        signal=struct;
        rx=struct;
        tx=struct;
        target=struct;
        noiseFloor; % noisefloor in dbm
        c = physconst('LightSpeed');
        k = physconst('Boltzmann');
        T=290;
    end
    
    methods
        function obj = c_Enviroment(envParam)
            obj.rx.coord=envParam.rx.coord;
            obj.tx.coord=envParam.tx.coord;
            obj.target.coord=envParam.target.coord;
            obj.target.cs=envParam.target.cs;
            obj.signal.cFq=envParam.signal.cFq;
            obj.signal.fs=envParam.signal.fs;
            obj.signal.pwr=envParam.signal.pwr;
            obj.signal.effBw=envParam.signal.effBw;
            obj.signal.noiseFig=1;
            
        end
        
        function out = propagate(obj,signal,dist)
            %function to simulate signal dampening due to dist
            %calc db attenuation
            att=((4*pi*dist*obj.signal.cFq)/obj.c)^2
            att=db2mag(10*log10(att))
            
            out=signal/att;
            
        end
        
        function out=powerDensityAtDist(obj,tPwr,dist)
           out=tPwr/(4*pi*dist^2) 
            
        end
        function out=propogateEcho(obj,signal,dist)
            lambda=obj.c/obj.signal.cFq;
            snr=(obj.signal.pwr*1*obj.target.cs*1*lambda^2)/((4*pi)^3*dist^4*obj.k*obj.T*obj.signal.effBw*obj.signal.noiseFig)
            mag2db(snr)
        end
        
        function out=powerDensityEcho(obj,tPwr,dist)
            out=(tPwr*1*obj.target.cs)/((4*pi)^2*dist^4)
        end
        
        function out = addNoiseFloor(obj,signal)
            %Todo add noise floor
        end
    end
end

