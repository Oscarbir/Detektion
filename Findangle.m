N = 1000;
s= randn(N,1)+1i*randn(N,1);
x1 = s + randn(N,1)+ 1i*randn(N,1);
plot([-2.2 2.2],[0 0],'k-'), hold on
plot([0 0],[-2.2 2.2],'k-')
axis([-2.2 2.2 -2.2 2.2]);
for j  = 1:20
    x2 = (s + randn(N,1)+1i*randn(N,1))*exp(j*1i*pi/10); % fasförkjutning
    r = xcorr(x1,x2,10)/length(x1);
    plot(real(r(11)),imag(r(11)),'r.','markersize',10)
    phi = rad2deg(angle(r(11)))
    pause(0.5)
end
 
 
 
 