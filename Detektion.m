clc, clear all
X=load('x.mat');
x=X.x;

%%
% x1=x(:,1);
% x2=x(:,2);
r1=xcorr(x(:,1),x(:,2),10000)/length(x(:,1));

 %subplot(1,4,1)
plot(abs(r1))

% [afmag,delay,doppler] = ambgfun(x1(1:200),x2(1:200),10e6,[100000 100000]);
% surf(delay*1e6,doppler/1e3,afmag,'LineStyle','none');
% axis tight; grid on; view([140,35]); colorbar;
% xlabel('Delay \tau (us)');ylabel('Doppler f_d (kHz)');
% title('Linear FM Pulse Waveform Ambiguity Function');

%%

y=fft(x(:,1));
t=linspace(1,1/10e6,10000);
plot(t,abs(y))


