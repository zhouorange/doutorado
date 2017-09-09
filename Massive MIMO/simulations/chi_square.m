clear all;clc

N = 1000;
M = 1000000;
a = (1/sqrt(2))*(randn(M,N) + 1i*randn(M,N));

chi = zeros(1,M);
for i=1:1:M
    chi(i) = sum(abs(a(i,:)).^2);
end

hist(chi,1000)

mean(chi)

var(chi)