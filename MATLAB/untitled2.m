A=zeros(1,10001)
for K=0:0.001:10

zB = 1 + 1/2*log(1 + pi*K^2/1024*(log(3125*pi*K^2/1024))^5);
A(int16(K*1000+1)) = 2/(zB*K)*(1-exp(-1/2*zB*K));

end

loglog(0:0.001:10,A)