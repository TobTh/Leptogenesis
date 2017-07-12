T = 15;
T1 = 40/sqrt(33/8 + 3*((3*pi)/2)^(2/3));
a = 33/8 + 3*((3*pi)/2)^(2/3);
b=7;
lambda=2;

x=1:0.1:50;

f1 = 1/2*a*(20^2 - T1^2)*x.^2 - 1/3*b*T*x.^3+1/4*lambda*x.^4; 
f2 = 1/2*a*(16^2 - T1^2)*x.^2-1/3*b*T*x.^3+1/4*lambda*x.^4; 
f3 = 1/2*a*(15^2 - T1^2)*x.^2-1/3*b*T*x.^3+1/4*lambda*x.^4; 
f4 = 1/2*a*(14.7^2 - T1^2)*x.^2-1/3*b*T*x.^3+1/4*lambda*x.^4;

plot(x,f1)
hold on
plot(x,f2)
plot(x,f3)
plot(x,f4)
hold off
xlabel('\phi')
ylabel('V_{eff}')