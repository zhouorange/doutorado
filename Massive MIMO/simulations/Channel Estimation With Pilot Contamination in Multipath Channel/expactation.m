clear all;close all;clc

theta_ik = 10;

N = 10000000;
y = zeros(1,N);
%y = 0;
for i=1:1:N
    
    z = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
    
    %y = y + (conj(z)*z)*(z*conj(z));
    
    y(i) = conj(z)*(z*conj(z));
    %y(i) = (z*conj(z))*(z*conj(z));
end
%y = y/N
%mean(y)

hist(real(y),10000)

figure
hist(imag(y),10000)