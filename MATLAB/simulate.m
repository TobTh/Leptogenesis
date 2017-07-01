clc
 
z = linspace(1,10^5);
y0 = [0 0 0];
y1 = [0 3/4 0];
x0 = [0 0];
x1 = [0,3/4];

 A1=zeros(1,91);
 A2=zeros(1,91);
 A3=zeros(1,91);
 A4=zeros(1,91);
 
 
for K=1:0.1:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K),z,y0);
    %[z2,N2] = ode15s(@(z,N) test3(z,N,K),z,y1);
    A1(1,int16(K*10-9)) = N1(end,1);
    %A2(1,int16(K*10-9)) = N2(end,1);
    [z3,N3] = ode15s(@(z,N) density_nonCorr(z,N,K),z,x0);
    [z4,N4] = ode15s(@(z,N) density_nonCorr(z,N,K),z,x1);
    A3(1,int16(K*10-9)) = N3(end,1);
    A4(1,int16(K*10-9)) = N4(end,1);
end


loglog((1:0.1:10),abs(A1.*(10^6*4/3)))
hold on
loglog((1:0.1:10),abs(A2.*(10^6*4/3)))
loglog((1:0.1:10),abs(A3.*(10^6*4/3)))
loglog((1:0.1:10),abs(A4.*(10^6*4/3)))
hold off