import Usefulfunctions.*

clear all
clc

Fs=100; % sample frekv
N=8192; % samples
Ts=1/100; %samle intervall
Tmax=(N-1)*Ts; 
t=0:Ts:Tmax; % tidsintervall

s=randn(N,1)+1i*rand(N,1);

bw=linspace(4000,1,6000)

Gain=30;
for n=bw:1:length(bw)
    %s=(randn(N,1)+1i*rand(N,1));
    s(n)=Gain*s(n);
end

s2=delay(s,1);
plot(t,10*log10(abs(s)))
