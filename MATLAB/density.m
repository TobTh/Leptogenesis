function ndot = density(z,N,K)        
    md = 1;%0^7;
    ms = 1;%.08*10^(-3);
    T = md/(K*ms)*1/z;
    M = 1;%0^10;
    cl = 1;
    cp = 2/3;
    l = 1;
    g1 = 1;
    g2 = 1;
    ht = 1;
    h11 = 1;%0;
    epsilon = 1;%0^(-6);
    Gamma0 = h11^2*M/(8*pi);
    H = 1/z^2;%*1.66*106.75*M^2/(1.22*10^28);
    nE = 1/pi^2*z^2*T^3*besselk(2,z);
    GammaBL =  3/pi^2*(cl + cp/2)*z^2*besselk(1,z)*Gamma0;
    a = 1 - l/z^2 - ht^2*(21/(2*(4*pi)^2) + 7*pi^2/(60*z^4))+...
        +(g1^2 + 3*g2^2)*(29/(8*(4*pi)^2) - pi^2/(80*z^4));
    b = -(ht^2*7*pi^2/45 + (g1^2 + 3*g2^2)*pi^2/60)/z^4;
    uE = 3/(2*pi^2)*T^3*besselk(3,z);
    
    ndot=[(epsilon*Gamma0*(N(2)-nE)-GammaBL*N(1)+...
               +epsilon*Gamma0*(N(3)-uE)-3*H*N(1))/(z*H);...
               (-a*Gamma0*(N(2)-nE)+...
               +(a-2*b)*Gamma0*(N(3)-uE)-3*H*N(2))/(z*H);...
               (-a*Gamma0*(N(3)-uE)-5*H*N(3))/(z*H)];
    
    
 
 
    