clc, clear all

%%T�nker att vi skapar en signal f�rst 
X=load('x.mat');
x=X.x;

r1=xcorr(x(:,1),x(:,2),3000)/length(x(:,1));
r2=xcorr(x(:,1),x(:,3),3000)/length(x(:,1));
r3=xcorr(x(:,1),x(:,4),3000)/length(x(:,1));
r4=xcorr(x(:,1),x(:,5),3000)/length(x(:,1));
 subplot(1,4,1)
plot(abs(r1))
 subplot(1,4,2)
 plot(abs(r2))
subplot(1,4,3)
plot(abs(r3))
subplot(1,4,4)
plot(abs(r4))

figure(2)

pspectrum(x(:,1))

%% genSignal
import Usefulfunctions.*
close,clear all,clc

sParam=struct;
sParam.centerFq=0;
sParam.bw=2e6;
sParam.gain=50;
sParam.fs=10e6;
sParam.N=10e3;
sParam.noise=0.1;

signalGen=c_SignalGen(sParam);
[fAxis,out]=signalGen.generateSignal;

pspectrum(out,sParam.fs)
%% genSignal and RFpropogate
clc, close, clear all
%init signal
sParam=struct;
sParam.centerFq=650e6;
sParam.bw=2e6;
sParam.power=100; %power in watt
sParam.fs=10e6; %samplefq
sParam.N=10e3;  %nr samples
sParam.noise=0; %noise 

%inspect signal..
signalGen=c_SignalGen(sParam);
[fAxis,out]=signalGen.generateSignal;
watt_signal=bandpower(out,sParam.fs,[-1e6 1e6])
figure(1)
pspectrum(out,sParam.fs)
hold on
% todo plotta amgfunc

%init enviroment
envParam=struct;
envParam.signal.cFq=sParam.centerFq;
envParam.signal.fs=sParam.fs;
envParam.signal.pwr=sParam.power;
envParam.signal.effBw=sParam.fs;
envParam.rx.coord=[0,0];
envParam.tx.coord=[5e3,0];
envParam.target.coord=[(2.5e3),(2.5e3)];
envParam.target.cs=40;

env=c_Enviroment(envParam);

out=env.propagate(out,5e3);
%env.propogateEcho(out,7070);

dPwr=env.powerDensityAtDist(sParam.power,5e3)
rPwr=env.powerDensityEcho(sParam.power,7070)

att=10*log10(rPwr/dPwr)


pspectrum(out,sParam.fs)
pspectrum(out*(db2mag(att)),sParam.fs)

