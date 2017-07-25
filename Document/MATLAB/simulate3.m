clc

a = 0.1646; %= 1/pi^2*besselk(2,1)
b = 1.0793; %= 3/(2*pi^2)*besselk(3,1);

z = linspace(1,20);
y1 = [0 a b]; 
x1 = [0,a];

A13=zeros(1,191);
A23=zeros(1,191);
A33=zeros(1,191);
A43=zeros(1,191);
A53=zeros(1,191);
A63=zeros(1,191);
A73=zeros(1,191);
A83=zeros(1,191);
A93=zeros(1,191);
A103=zeros(1,191);

cl = 1;
cp = 2/3;

 
for K=1:0.1:20
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A13(1,int16(K*10-9)) = N1(end,1);
    A23(1,int16(K*10-9)) = N2(end,1);
end


cl = 1;
cp = 14/23;

 
for K=1:0.1:20
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A33(1,int16(K*10-9)) = N1(end,1);
    A43(1,int16(K*10-9)) = N2(end,1);
end


cl = 3/4;
cp = 1/2;
 
for K=1:0.1:20
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A53(1,int16(K*10-9)) = N1(end,1);
    A63(1,int16(K*10-9)) = N2(end,1);
end


cl = 78/115;
cp = 56/115;

for K=1:0.1:20
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A73(1,int16(K*10-9)) = N1(end,1);
    A83(1,int16(K*10-9)) = N2(end,1);
end


cl = 344/537;
cp = 52/179;

for K=1:0.1:20
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A93(1,int16(K*10-9)) = N1(end,1);
    A103(1,int16(K*10-9)) = N2(end,1);
end

quinticMA1 = sgolayfilt((A13-A23)./A13, 10, 101);
quinticMA2 = sgolayfilt((A33-A43)./A33, 10, 101);
quinticMA3 = sgolayfilt((A53-A63)./A53, 10, 101);
quinticMA4 = sgolayfilt((A73-A83)./A73, 10, 101);
quinticMA5 = sgolayfilt((A93-A103)./A93, 10, 101);

figure(2)
 semilogx((1:0.1:20),quinticMA1)
 hold on
 semilogx((1:0.1:20),quinticMA2) 
 semilogx((1:0.1:20),quinticMA3) 
 semilogx((1:0.1:20),quinticMA4)
 semilogx((1:0.1:20),quinticMA5)
 hold off
 xlabel('K')
 ylabel('(\kappa-\kappa_{NR})/\kappa')


function ndot = density_Rel(z,N,K,cl,cp)        
    nE = 1/pi^2*z^2*besselk(2,z);
    uE = 3/(2*pi^2)*besselk(3,z)*z^3;
   
    ndot=[5.4737*z*K*((N(2)-nE)-1/z^2*(N(3)-uE))-K*z^3*besselk(1,z)*3/pi^2*(cl+cp/2)*N(1);...
        -K*z*(N(2)-nE)+K/z*(N(3)-uE);...
        -K*(N(3)-uE)*z];
end

function ndot = density_nonCorr(z,N,K,cl,cp)               
    nE = 1/pi^2*z^2*besselk(2,z);
    
    ndot=[z*5.4737*K*(N(2)-nE)-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*(N(2)-nE)*z];
end



