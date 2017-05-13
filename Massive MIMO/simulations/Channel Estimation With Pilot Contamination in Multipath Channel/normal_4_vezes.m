clear all;clc;close all;

N = 1000000;

res = zeros(1,N);
for i=1:1:N
    
    %x = randn(1,1);
    
    x = (1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));
    
    res(i) = (x'*x)*(x'*x);
    
end

mean(res)

var(res)