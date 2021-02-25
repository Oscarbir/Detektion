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
    end
    
    methods
        function obj = c_Enviroment(envParam)
            obj.rx.coord=envParam.rx.coord;
            obj.tx.coord=envParam.tx.coord;
            obj.target.coord=envParam.target.coord;
            obj.signal.cFq=envParam.signal.cFq;
            obj.signal.fs=envParam.signal.fs;
            
        end
        
        function out = propagate(obj,signal,dist)
            %function to simulate signal dampening due to dist
            %calc db attenuation
            att=((4*pi*dist*obj.signal.cFq)/obj.c)^2
            att=db2mag(10*log10(att))
            
            out=signal/att;
            
        end
    end
end

