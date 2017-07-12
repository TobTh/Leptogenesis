clc
 
a = 9.7046/90; %=45/4*zeta(5)/zeta(3)

z = linspace(1,10^6);
y1 = [0 3/4 a];
x1 = [0,3/4];

A13=zeros(1,901);
A23=zeros(1,901);
A33=zeros(1,901);
A43=zeros(1,901);
A53=zeros(1,901);
A63=zeros(1,901);
A73=zeros(1,901);
A83=zeros(1,901);
A93=zeros(1,901);
A103=zeros(1,901);

cl = 1;
cp = 2/3;

 
for K=1:0.01:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A13(1,int16(K*100-99)) = N1(end,1);
    A23(1,int16(K*100-99)) = N2(end,1);
end


cl = 1;
cp = 14/23;

 
for K=1:0.01:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A33(1,int16(K*100-99)) = N1(end,1);
    A43(1,int16(K*100-99)) = N2(end,1);
end


cl = 3/4;
cp = 1/2;
 
for K=1:0.01:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A53(1,int16(K*100-99)) = N1(end,1);
    A63(1,int16(K*100-99)) = N2(end,1);
end


cl = 78/115;
cp = 56/115;

for K=1:0.01:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A73(1,int16(K*100-99)) = N1(end,1);
    A83(1,int16(K*100-99)) = N2(end,1);
end


cl = 344/537;
cp = 52/179;

for K=1:0.01:10
    
    [z1,N1] = ode15s(@(z,N) test3(z,N,K,cl,cp),z,y1);
    [z2,N2] = ode15s(@(z,N) density_nonCorr(z,N,K,cl,cp),z,x1);
    A93(1,int16(K*100-99)) = N1(end,1);
    A103(1,int16(K*100-99)) = N2(end,1);
end

loglog((1:0.01:10),abs((A13-A23)./A13))
hold on
loglog((1:0.01:10),abs((A33-A43)./A33))
loglog((1:0.01:10),abs((A53-A63)./A53))
loglog((1:0.01:10),abs((A63-A83)./A73))
loglog((1:0.01:10),abs((A73-A103)./A93))
hold off

% loglog((1:0.01:10),abs(A13.*(10^6*4/3)))
% hold on
% loglog((1:0.01:10),abs(A23.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A33.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A43.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A53.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A63.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A73.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A83.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A93.*(10^6*4/3)))
% loglog((1:0.01:10),abs(A103.*(10^6*4/3)))
% hold off