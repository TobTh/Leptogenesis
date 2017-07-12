clc

cl = 1;
cp = 0;

z = linspace(1,10^3);
x0 = [0 0];
x1 = [0,3/4];

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


semilogx((1:0.5:100),abs(A32./A12))
hold on
semilogx((1:0.5:100),abs(A42./A22))
hold off
xlabel('K')
ylabel('\kappa/\kappa_{classical}')


% loglog((1:0.01:10),abs(A1.*(10^6*4/3)))
% hold on
% loglog((1:0.01:10),abs(A2.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A3.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A4.*(10^6*4/3)))
% hold off