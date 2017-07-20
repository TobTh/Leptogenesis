clc
 

cl = 344/537;
cp = 52/179;

a = 0.1646; %= 1/pi^2*besselk(2,1)
b = 1.0793; %= 3/(2*pi^2)*besselk(3,1);

z = linspace(1,20);
y0 = [0 0 0];
y1 = [0 a b]; 
x0 = [0 0];
x1 = [0,a];

 A11=zeros(1,37);
 A21=zeros(1,37);
 A31=zeros(1,37);
 A41=zeros(1,37);
 
 
for K=1:0.25:10
    
    [z1,N1] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y0);
    [z2,N2] = ode15s(@(z,N) density_Rel(z,N,K,cl,cp),z,y1);
    A11(1,int16(K*4-3)) = N1(end,1);
    A21(1,int16(K*4-3)) = N2(end,1);
    [z3,N3] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x0);
    [z4,N4] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A31(1,int16(K*4-3)) = N3(end,1);
    A41(1,int16(K*4-3)) = N4(end,1);
end

quinticMA1 = sgolayfilt(A11, 10, 31);
quinticMA2 = sgolayfilt(A21, 10, 31);
quinticMA3 = sgolayfilt(A31, 10, 31);
quinticMA4 = sgolayfilt(A41, 10, 31);

figure(1)
 loglog((1:0.25:10),quinticMA1)
 hold on
 loglog((1:0.25:10),quinticMA2)
 loglog((1:0.25:10),quinticMA3)
 loglog((1:0.25:10),quinticMA4)
 hold off
 xlabel('K')
 ylabel('\kappa')

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