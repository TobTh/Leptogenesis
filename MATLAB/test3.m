function ndot = test3(z,N,K,cl,cp)        
    epsilon = 10^(-6);
    nE = 3/8*z^2*besselk(2,z);
    uE = 0.8319*besselk(3,z)*z;
   
    ndot=[z*(epsilon*K*(N(2)-nE))+epsilon*K*z*(N(3)-uE)-3/pi^2*(cl + cp/2)*z^3*besselk(1,z)*K*N(1);...
        -K*(N(2)-nE)*z+epsilon*K*z*(N(3)-uE);...
        -K*(N(3)-uE)*z];