function ndot = density_nonCorr_classic(z,N,K,cl,cp)               
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    
    ndot=[z*(epsilon*K*(N(2)-nE))-1/4*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*(N(2)-nE)*z];
    
 
 
    