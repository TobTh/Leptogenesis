function ndot = density_nonRel(z,N,K)        
    cl = 344/537;
    cp = 52/179;
    l = 0.5;
    g1 = 1;
    g2 = 1;
    ht = 1;
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    a = 1 - l/z^2 - ht^2*(21/(2*(4*pi)^2) + 7*pi^2/(60*z^4))+...
        +(g1^2 + 3*g2^2)*(29/(8*(4*pi)^2) - pi^2/(80*z^4));
    
    ndot=[z*(epsilon*a*K*(N(2)-nE))-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
         -K*a*(N(2)-nE)*z];
    
    
 
 
    