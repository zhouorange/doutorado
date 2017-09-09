clear all;clc

N = 1000000;

M = 10;

ZZ = zeros(1,N);
for k=1:1:N
    z = 0;
    for i=1:1:M
        x = (1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));
        y = (1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));
        z = z + x'*y;
    end
    ZZ(k) = z;
end

%hist(ZZ,1000)

mean(ZZ)

var(ZZ)