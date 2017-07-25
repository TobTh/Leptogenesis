z=1:0.1:10;

f1=z.^2;
f2=3/pi^2*(1 + 1/3)*z.^4.*besselk(1, z); 
f3= 3/pi^2*(344/537 + 26/179)*z.^4.*besselk(1, z);

figure(5)
loglog(z,f1)
hold on
loglog(z,f2)
loglog(z,f3)
hold off
xlabel('z')
ylabel('K^{-1} \Gamma/H')