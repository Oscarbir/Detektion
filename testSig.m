clear all 
close
clf
clc

f=1e6;
fs=10e6;
dfs=1/fs;
N=10000;

k=0:1:N;


t=0:dfs:N*dfs;

wk=2*pi*k*fs/N; % rad-vektor
w=2*pi*f;

x=sin(wk*t)

plot(wk,fftshift(x))
%% sa
close 
clear all
clc

fc =40;
fs =9000;

wc = 2*pi*fc;
ws = 2*pi*fs;

t = 0:(1/fs):2;
x = sin(wc*t);

figure()
plot(t,x)
n =length(t);
f = (0:n-1)*(fs/n);
N =1000;
F = fft(x,N)/n;
f = (-N/2:N/2-1)*(fs/N);
F0 =fftshift(F);

figure()
stem([-40, 40], [pi, pi]) % plot the spectrum of cos(w_c*t)
hold on
plot(f, abs(F0))
legend('spectrum of original signal','spec of sampled sig')
xlabel('f （HZ（')
ylabel('X(jw), (Hz)')
title('spectrum of sampled signal')
%% test
clear all 
close 
clc
import Usefulfunctions.*
clf

fs=1e6;

N=10000;
s= randn(N,1)+1i*randn(N,1);

Ns = 10;

x1 = s + randn(N,1)+ 1i*randn(N,1);
x2 = delay(s + randn(N,1)+1i*randn(N,1),-5);

r = xcorr(x1,x2,Ns)/length(x1);
plot(abs(r)), hold on 

x2 = delayC(x1,x2,'time');
r = xcorr(x1,x2,Ns)/length(x1);
plot(abs(r))
%% noise signal broadband
clear all 
close 
clc
import Usefulfunctions.*
clf

N=10e6;
fs=10e6;
bw=2e6;

df=N/fs

Nbw=bw/df


gain=30;

dt=1/fs;

s1=randn(N,1)+ 1i*randn(N,1);
s2=randn(Nbw,1)+ 1i*randn(Nbw,1);

vec= [1 ;2 ;3 ;4]

z=[zeros(N/2-Nbw/2,1); s2+gain ; zeros(N/2-Nbw/2,1,1)];

f = (-N/2:N/2-1)*(fs/N);

s3=z+s1;

df=fs/N
figure(1)

figure(2)
plot(f,abs(s3))
hold on
%plot(f,abs(z))
% 
% 


%% ss
[M,I] = max(r);

tau = I-Ns; 
x2 = delay(x2,tau);
r = xcorr(x1,x2,Ns)/length(x1);
plot(abs(r))

%% ads

x2 = delay(s + randn(N,1)+1i*randn(N,1))*exp(1i*pi/4); % fasförkjutning

r = xcorr(x1,x2,10)/length(x1);

plot([real(r),imag(r)]); % plottar real- och im-del för sig

%% frekvens-offset
import Usefulfunctions.*
clf

N=10000;
s= randn(N,1)+1i*randn(N,1);
t = (1:length(s))';
x1 = s + randn(N,1)+ 1i*randn(N,1);
x2 = 2*delay((s + randn(N,1)+1i*randn(N,1)),-5).*exp(1i*2*pi*0.0001*t); % frekvensoffset
%x2 = x2.*exp(-1i*2*pi*0.0001*t); 

% x1 = x1/std(x1);
% x2 = x2/std(x2);

% x2 = delayC(x1,x2,'time');
% x2 = delayC(x1,x2,'freq');
r = xcorr(x1,x2,10)/length(x1);
plot(abs(r))

%% sad
x2 = delayC(x1,x2,'freq');


%% asa
X1 = fft(x1);
X2 = fft(x2);

Ns = 10;
r = xcorr(X1,X2,Ns)/length(X1);

X2 = delayC(X1,X2);
r = xcorr(X1,X2,Ns)/length(X1);

x2 = ifft(X2);
r = xcorr(x1,x2)/length(x1);
plot(abs(r))

%% fsf
% skapa nu signal som kompenserar

x2_test = x2.*exp(-1i*2*pi*0.0001*t); 
% testa olika frekvenser till det matchar (tills max är störst) 

r = xcorr(X1,X2)/length(x1);

plot(abs(r))

%% sff

% xcorr för FFT

% Maximera funktionen 

% fMinSearch

% Hitta integer-offset för frekvens och tidsdomän

% IQ-mixer

% Falska mål pga IQ-mixern



