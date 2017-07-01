function ndot = test3(z,N,K)        
    M = 1;
    cl = 344/537;
    cp = 52/179;
    l = 0.5;
    g1 = 1;
    g2 = 1;
    ht = 1;
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    uE = 1/(5*zeta(3))*zeta(7/2)/sqrt(pi)*2^3*6*besselk(3,z)*M^3;
   
    ndot=[z*(-epsilon*K*(N(2)-nE))+epsilon*K*z*(N(3)-uE)-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
        -K*(N(2)-nE)*z+epsilon*K*z*(N(3)-uE);...
        -K*(N(3)-uE)*z];