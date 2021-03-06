clc, clear all
X=load('x.mat');
x=X.x;

%%
% x1=x(:,1);
% x2=x(:,2);

 %subplot(1,4,1)

r1=xcorr(x(:,1),x(:,2),10000)/length(x(:,1));
r2=xcorr(x(:,1),x(:,3),100)/length(x(:,1));
r3=xcorr(x(:,1),x(:,4),100)/length(x(:,1));
r4=xcorr(x(:,1),x(:,5),100)/length(x(:,1));
phi = rad2deg(angle(r1(9033)))  

%  subplot(1,4,1)
plot(abs(r1))
hold on 
plot(
%  subplot(1,4,2)
% plot([real(r2),imag(r2)])
%  subplot(1,4,3)
% plot([real(r3),imag(r3)])
%  subplot(1,4,4)
% plot([real(r4),imag(r4)])
% plot(abs(r1))

% [afmag,delay,doppler] = ambgfun(x1(1:200),x2(1:200),10e6,[100000 100000]);
% surf(delay*1e6,doppler/1e3,afmag,'LineStyle','none');
% axis tight; grid on; view([140,35]); colorbar;
% xlabel('Delay \tau (us)');ylabel('Doppler f_d (kHz)');
% title('Linear FM Pulse Waveform Ambiguity Function');


%%

y=fft(x(:,1));
t=linspace(1,1/10e6,10000);
plot(t,abs(y))



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
sParam=struct;
sParam.centerFq=0;
sParam.bw=2e6;
sParam.power=100; %power in watt
sParam.fs=10e6; %samplefq
sParam.N=10e3;  %nr samples
sParam.noise=0.1;


signalGen=c_SignalGen(sParam);
[fAxis,out]=signalGen.generateSignal;

watt_signal=rms(out)^2

pspectrum(out,sParam.fs)
