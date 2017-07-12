clc
 
%a = 9.7046; %=45/4*zeta(5)/zeta(3)
%a = 17.7228; %=3/zeta(3)*besselk(3,1)
a=0.8319;

cl = 344/537;
cp = 52/179;

z = linspace(1,10^6);
y0 = [0 0 0];
y1 = [0 3/4 a];
x0 = [0 0];
x1 = [0,3/4];

 A11=zeros(1,91);
 A21=zeros(1,91);
 A31=zeros(1,91);
 A41=zeros(1,91);
 
 
for K=1:0.1:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y0);
    [z2,N2] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    A11(1,int16(K*10-9)) = N1(end,1);
    A21(1,int16(K*10-9)) = N2(end,1);
    [z3,N3] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x0);
    [z4,N4] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A31(1,int16(K*10-9)) = N3(end,1);
    A41(1,int16(K*10-9)) = N4(end,1);
end


loglog((1:0.1:10),A11.*(10^6*4/3))
hold on
loglog((1:0.1:10),A21.*(10^6*4/3))
loglog((1:0.1:10),A31.*(10^6*4/3))
loglog((1:0.1:10),A41.*(10^6*4/3))
hold off