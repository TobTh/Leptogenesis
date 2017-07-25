clc
 

cl = 344/537;
cp = 52/179;

a = 0.1646; %= 1/pi^2*besselk(2,1)
b = 1.0793; %= 3/(2*pi^2)*besselk(3,1);

z = linspace(1,20);     % Range of z over which rate equations are solved
y0 = [0 0 0];           % Zero initial conditions for X_(B-L), X_N, X_u
y1 = [0 a b];           % Thermal intial conditions for X_(B-L), X_N, X_u
x0 = [0 0];             % Zero initial conditions for X_(B-L), X_N
x1 = [0,a];             % Thermal intial conditions for X_(B-L), X_N

 A11=zeros(1,37);       % Matrices of zeroes for later storage of data
 A21=zeros(1,37);
 A31=zeros(1,37);
 A41=zeros(1,37);
 
 
% The rate equations are solved for each value of K in a range from 1 to 10
% in steps of 0.25. 
% The last calculated data point for each K gets stored in the matrices
% defined above
 
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

% Obtained results get smoothed out in order to suppress noise
 
quinticMA1 = sgolayfilt(A11, 10, 31);
quinticMA2 = sgolayfilt(A21, 10, 31);
quinticMA3 = sgolayfilt(A31, 10, 31);
quinticMA4 = sgolayfilt(A41, 10, 31);

% Smoothed results are plotted against K

figure(1)
 loglog((1:0.25:10),quinticMA1)
 hold on
 loglog((1:0.25:10),quinticMA2)
 loglog((1:0.25:10),quinticMA3)
 loglog((1:0.25:10),quinticMA4)
 hold off
 xlabel('K')
 ylabel('\kappa')


 
 % This function defines the rate equation with relativistic corrections
 % included
 
 function ndot = density_Rel(z,N,K,cl,cp)        
    nE = 1/pi^2*z^2*besselk(2,z); % X_N in equilibrium 
    uE = 3/(2*pi^2)*besselk(3,z)*z^3; % X_u in equilibrium 
    
    ndot=[5.4737*z*K*((N(2)-nE)-1/z^2*(N(3)-uE))-K*z^3*besselk(1,z)...
        *3/pi^2*(cl+cp/2)*N(1);...
        -K*z*(N(2)-nE)+K/z*(N(3)-uE);...
        -K*(N(3)-uE)*z];
 end
 

 
% This function defines the rate equation without relativistic corrections
% included

function ndot = density_nonCorr(z,N,K,cl,cp)               
    nE = 1/pi^2*z^2*besselk(2,z); % X_N in equilibrium 
    
    ndot=[z*5.4737*K*(N(2)-nE)-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*...
        K*N(1);...
         -K*(N(2)-nE)*z];
end