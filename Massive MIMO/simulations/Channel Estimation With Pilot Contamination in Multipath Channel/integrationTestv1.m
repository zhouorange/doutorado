clear all;close all;clc

rng(1)

M = 10;            % Number of antennas.
K = 2;              % Number of single-antenna users.
P = 4;              % Channel Length (Finite Impulse Response - FIR).
L = 2;              % Number of cells.

N = getPilotLength(K,P);  % Pilot length is set according to K and P.
q = 100000000;            % Uplink pilot power.

beta = abs(randn(L,K));
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);
epsilon11 = (beta_sum + 1/(q*N));

kik = calculateKik(beta111, epsilon11);

fun = @(w,t) (((kik.^2).*(1-t)+kik.*w.*sqrt(t.*(1-t)))./((kik.^2).*(1-t)+2.*kik.*w.*sqrt(t.*(1-t))+t)).*((gamma(2.*M)./(gamma(M).^2)).*((t.*(1-t)).^(M-1))).*((M./pi).*((1-(w.^2)).^(M-0.5)));

q = integral2(fun,-1,1,0,1,'Method','iterated','AbsTol',0,'RelTol',1e-5);