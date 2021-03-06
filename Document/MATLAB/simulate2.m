clc

cl = 1;
cp = 0;

a = 0.1646; %= 1/pi^2*besselk(2,1)

z = linspace(1,20);
x0 = [0 0]; 
x1 = [0,a];


A12=zeros(1,199);
A22=zeros(1,199);
A32=zeros(1,199);
A42=zeros(1,199);
 
 
for K=1:0.5:100
    
    [z1,N1] = ode15s(@(z,N) density_nonCorr_classic(z,N,K,cl,cp),z,x0);
    [z2,N2] = ode15s(@(z,N) density_nonCorr_classic(z,N,K,cl,cp),z,x1);
    A12(1,int16(K*2-1)) = N1(end,1);
    A22(1,int16(K*2-1)) = N2(end,1);
    [z3,N3] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x0);
    [z4,N4] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A32(1,int16(K*2-1)) = N3(end,1);
    A42(1,int16(K*2-1)) = N4(end,1);
end

figure(3)
 quinticMA1 = sgolayfilt(A32./A12, 10, 31);
quinticMA2 = sgolayfilt(A42./A22, 10, 31);

 semilogx((1:0.5:100),quinticMA1)
 hold on
 semilogx((1:0.5:100),quinticMA2)
 hold off
 xlabel('K')
 ylabel('\kappa/\kappa_{classical}')

function ndot = density_nonCorr(z,N,K,cl,cp)               
    nE = 1/pi^2*z^2*besselk(2,z);
    
    ndot=[z*5.4737*K*(N(2)-nE)-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*(N(2)-nE)*z];
end

function ndot = density_nonCorr_classic(z,N,K,cl,cp)               
    nE = 1/pi^2*z^2*besselk(2,z);
    
    ndot=[z*5.4737*K*(N(2)-nE)-1/4*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*(N(2)-nE)*z];
end