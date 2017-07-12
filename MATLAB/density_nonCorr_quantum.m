function ndot = density_nonCorr_quantum(z,N,K,cl,cp)        
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    %fun = sqrt(x.^2-1)/(exp(x.*z-1));
    %I = 4*pi*integral(fun,1,10^1000);
    
    ndot=[z*(epsilon*K*(N(2)-nE))-3/pi^2*z^4*exp(z/2)/(exp(z/2)-1)*(cp/2+cl)...
        *integral(@(x) sqrt(x.^2-1)./(exp(x*z)-1),1,10^1000)*K*N(1);...
         -K*(N(2)-nE)*z];
    
    
 
 
    