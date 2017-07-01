function ndot = density_nonCorr(z,N,K)               
    cl = 344/537;
    cp = 52/179;
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    
    ndot=[z*(epsilon*K*(N(2)-nE))-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*(N(2)-nE)*z];
    
 
 
    