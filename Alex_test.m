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
