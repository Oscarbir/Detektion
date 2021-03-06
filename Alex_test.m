clc, clear all,clf

Tx=[-2320 0]; %%Sändare
Rx=[2320 0];     %%Mottagare

%Targ = rand(1,2);   %% Target på en godtycklig koordinat
%Targ(1)= Targ(1)*10000;
%Targ(2)= Targ(2)*10000;
%%Geometri
%R1 = sqrt((Targ(1)-Tx(1)).^2+(Targ(2)-Tx(2)).^2);
%R2 = sqrt((Targ(1)-Rx(1)).^2+(Targ(2)-Rx(2)).^2);
[x1,y1] = pol2cart(deg2rad(31),15000)

plot([2320,x1],[0,y1],'k-');
hold on
Rb = sqrt((Tx(1)-Rx(1)).^2+(Tx(2)-Rx(2)).^2);

%BR = R1 + R2 - Rb; %Bistatisk range

Re=29000-4800;

a = 0.5*abs(Re);
b = sqrt((Re).^2-Rb.^2)./2;


radie=10000;
t=-pi:0.01:pi;
x=a*cos(t);
y=b*sin(t);
plot(x,y)
hold on 
plot(2320+radie*cos(t),radie*sin(t))
plot(Tx(1),Tx(2),'.','markersize',20)
plot(Rx(1),Rx(2),'.','markersize',20)


%plot(Targ(1),Targ(2),'.','markersize',20)
hold off


%% elippsss

clear all
clf

Rl=4.7e3;

signaldelay=5e-6;
distdelay=signaldelay*physconst('lightspeed');

f1x = -Rl/2;
f1y = 0
f2x = Rl/2;
f2y = 0
%bistaticrange = Rt+Rr-Rl, ellipsesum = Rt+Rr. 

ellipsesum = distdelay+Rl % constant sum 

a=ellipsesum/2; %major axis
x0 = (f1x+f2x)/2; %center coord
y0 = (f1y+f2y)/2;
f = sqrt((f1x-x0)^2+(f1y-y0)^2) %dist from focal to center
b = sqrt(a^2-f^2); %minor axis


t = linspace(0, 2*pi, 1000)
x = a*cos(t); 
y = b*sin(t);

plot(x,y)
hold on
plot(f1x,f1y,'x')
plot(f1x,f2y,'o')
axis equal


for n=1:length(t)
    dx=x(n)-f1x;
    dy=y(n)-f1y;
    rad=atan2(dx,dy);
    deg(n)=rad2deg(rad);
    deg(n) = mod((deg(n) + 360), 360);
end

bearing=0;

[diff,n]=min(abs(deg-bearing))

Rr=sqrt((x(n)-f1x)^2+(y(n)-f1y)^2)
Rt=sqrt((x(n)-f2x)^2+(y(n)-f2y)^2)

(Rt+Rr-Rl)/3e8

figure(1)
plot([f1x x(n)],[f1y y(n)])
plot([f2x x(n)],[f2y y(n)])
