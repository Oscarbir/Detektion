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

