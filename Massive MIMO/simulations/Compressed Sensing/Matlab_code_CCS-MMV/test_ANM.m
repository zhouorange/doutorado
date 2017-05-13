% test of ANM

clear all
close all
clc

N = 20;
M = 10;
K = 3;
L = 5;
p = randperm(N);
Omega = sort(p(1:M))';

% true frequency and power
f0 = [.1; .2; .5];
A0 = exp(1i*2*pi*kron((0:N-1)',f0'));
S0 = (randn(K,L) + 1i*randn(K,L)) / sqrt(2);

% noiseless data
Y0 = A0 * S0;

% noisy samples
sigma2 = .01;
E = sqrt(sigma2/2) * (randn(M,L)+1i*randn(M,L));
YonOmega = Y0(Omega,:) + E;
eta = norm(E,'fro');

% CCS
tic
[Y,u,freq,amp] = ANM_sdpt3(YonOmega, Omega, N, eta);
toc

% MSE
if length(freq) < K
    MSE = inf;
else
    MSE = sqrt(norm(f0 - sort(freq(1:K)))^2 / K);
end

% plot result
figure(100), plot(f0, sqrt(mean(abs(S0).^2,2)), 'ko'), hold on;
plot(freq, amp, 'rx', 'MarkerSize', 5);
xlabel('Frequency');
ylabel('Power');
legend('Truth', 'Estimate');
title(strcat('Frequency MSE = ', num2str(MSE)));