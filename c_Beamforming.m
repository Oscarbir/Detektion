clc, clear all, close 
addpath /home/marc/Skola/Kandidat/antennaFunctions; 

%%simple antenna diagram xy

d=0.25;
a=[1,1];
[g, phi] = gain1d(d, a, 400);

dbz(phi, g, 30, 20)

%% asd
close, clear
d=0.5; N=8;
a = uniform(d, 90, N);
[g, phi] = gain1d(d, a, 400);
A = sqrt(g);
psi = 2*pi*d*cos(phi);
plot(psi/pi, A);
figure(2);
dbz(phi, g, 45, 20);

%% steeriing
d=0.5; N=11; ph0=60;
a = uniform(d, ph0, N);
[g, phi] = gain1/d(d, a, 400);
psi = 2*pi*d*cos(phi);
figure; plot(psi/pi, sqrt(g));
figure; dbz(phi, g, 30, 20);

%% test 
clc,clf
X=load('x.mat');
x=X.x;
figure(1)

pspectrum(sum(x,5))
a=[1 1 1 1 1];

lambda = physconst('LightSpeed')/(650e6)
d=0.5*lambda;
kd=2*pi*d/lambda;

%x(:,1)=x(:,1)*exp(1i*kd*cos(pi/2));
x(:,2)=a(2)*x(:,2).*exp(1i*kd*2*cos(0));
x(:,3)=a(3)*x(:,3).*exp(1i*kd*3*cos(pi/3)+pi/2);
x(:,4)=a(4)*x(:,4).*exp(1i*kd*4*cos(pi/3)+pi/2);
x(:,5)=a(5)*x(:,5).*exp(1i*kd*5*cos(pi/3)+pi/2);

Xs=x(:,1)+x(:,2)+x(:,3)+x(:,4)+x(:,5);
hold on
pspectrum(Xs)