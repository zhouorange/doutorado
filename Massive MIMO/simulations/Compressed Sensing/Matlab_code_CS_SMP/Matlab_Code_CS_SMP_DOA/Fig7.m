% This code is used to generate Fig. 7 of our paper:
% Z. Yang, C. Zhang, and L. Xie, "Robustly stable signal recovery 
% in compressed sensing with structured matrix perturbation," 
% IEEE Transactions on Signal Processing, vol. 60, no. 9, pp. 4658--4671, 2012.


close all;
clear all
clc

M = 30;
K = 2;

% true DOA
theta = [3/90; 13/90] + 2 * (rand(2,1)-.5) / 90;

% true Phi
Phi = zeros(M,K);
for m = 1:M
    for k = 1:K
        Phi(m,k) = exp(1i * pi * (m-(M+1)/2) * theta(k));
    end
end

% true signal
phase = rand(2,1);
X = exp(1i*2*pi*phase);

% observation
Y = Phi * X;



%%%%%%%%%%%%%%%% CS-SMP for off-grid DOA estimation %%%%%%%%%%%%
N = 90;
resolution = 2 / N;
grid = (-1+resolution/2:resolution:1-resolution/2)';

% uniform linear array (ULA), with the origin at the middle
A = zeros(M,N);
B = zeros(M,N);
for m = 1:M
    for n = 1:N
        temp = exp(1i * pi * (m-(M+1)/2) * grid(n));
        A(m,n) = temp;
        B(m,n) = 1i * pi * (m-(M+1)/2) * temp;
    end
end

beta = zeros(N,1); 
beta(N/2 + (round(theta/resolution+.5))) = theta - grid(N/2 + (round(theta/resolution+.5)));
X0 = zeros(N,1); X0(N/2 + (round(theta/resolution+.5))) = X; 

epsilon = pi^2 * sqrt(M/15 * (3*M^4-10*M^2+7)) / (4 * N^2);

maxiter = 200;
tol = 1e-6;
tstart = tic;
res_CS_SMP = AA_P_BPDN_DOA(Y, A, B, resolution/2, epsilon, maxiter, tol);
time = toc(tstart);

xp_rec1 = grid + res_CS_SMP.beta;


kappa = pi/2 * sqrt((M^2-1)/3);




%%%%%%%%%%%%% standard CS for on-grid DOA estimation %%%%%%%%%%%%
N = 360;

resolution = 2 / N;
grid = (-1+resolution/2:resolution:1-resolution/2)';
X1 = zeros(N,1); X1(N/2 + (round(theta/resolution+.5))) = X; 

% uniform linear array (ULA), with the origin at the middle
A = zeros(M,N);
for m = 1:M
    for n = 1:N
        temp = exp(1i * pi * (m-(M+1)/2) * grid(n));
        A(m,n) = temp;
    end
end

epsilon = pi / N * sqrt(M*(M^2-1)/3);

x = zeros(N,1);
cvx_begin
  variable x(N) complex;
  minimize ( norm(x,1) )
    subject to
      norm(Y - A * x) <= epsilon;
cvx_end


support = find(beta ~= 0);

figure(5),
subplot(2,2,1), stem((1:90)', abs(X0), 'ko', 'MarkerSize', 5);hold on;
stem((1:90)', abs(res_CS_SMP.x),'r*', 'MarkerSize', 5);
xlabel('Grid index'); ylabel('Signal amplitude');
axis([0,90,0,1]);
% legend('Original', 'Recovered');
subplot(2,2,2), stem((1:90)', beta * kappa, 'ko', 'MarkerSize', 5); hold on;
stem(support, res_CS_SMP.beta(support)*kappa,'r*', 'MarkerSize', 5);
xlabel('Grid index'); ylabel('\beta');
axis([0,90,-kappa/90,kappa/90]);
% legend('Original', 'Recovered');
subplot(2,2,3), stem((1:N)', abs(X1), 'ko', 'MarkerSize', 5);hold on;
stem((1:N)', abs(x),'b^', 'MarkerSize', 5);
xlabel('Grid index'); ylabel('Signal amplitude');
axis([0,360,0,1]);
% legend('Original', 'Recovered');
subplot(2,2,4), stem(theta,ones(K,1),'ko', 'MarkerSize', 5); hold on;
stem(xp_rec1,abs(res_CS_SMP.x),'r*', 'MarkerSize', 5);hold on;
stem(grid,abs(x),'b^', 'MarkerSize', 5);
axis([theta(1)-1.5/90,theta(1)+1.5/90,0,1]);
xlabel('Cosine DOA (near source 1)'); ylabel('Signal amplitude');
legend('Original', 'Recovered by SP-CS', 'Recovered by SCS');