clear all;close all;clc

M = 10;

n = M;
u = -0.4:0.01:0.4;
fu = zeros(1,length(u));
for u_idx=1:1:length(u)
    bessel1 = besseli(n,(1./(4*(u(u_idx).^2))));
    bessel2 = besseli(n+1,(1./(4*(u(u_idx).^2))));
    fu(u_idx) = exp(-1./(4*(u(u_idx).^2)))*( ( (8*n*u(u_idx).^2 - 1)*bessel1 + bessel2 ) / (4*(abs(u(u_idx)).^3)) );
end

plot(u,fu)
