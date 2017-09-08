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

theta_ik = calculateTheta_ik(beta111, epsilon11, M);

