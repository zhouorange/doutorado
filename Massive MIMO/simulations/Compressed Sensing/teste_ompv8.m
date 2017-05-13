clear all;close all;clc;

N = 2048;               % number of OFDM subcarriers.
NCP = 128;              % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
Fs = 30.72e6;
T = 1/Fs;

%% ********* Pilot sequence. *********
% Comment: The vector P addresses pilot sequences.
M = 32;                                         % Number of pilots. Length of the training sequence.
% Position of the pilot in the frequency domain.
%ppos = [482 526 579 669 690 738 932 1007 1018 1074 1075 1401 1543 1828 1878 1895]; % �tima distribui��o de frequ�ncia (subcarriers) para 16 pilotos.
%ppos = [42 92 121 185 375 531 892 1058 1224 1384 1466 1490 1544 1753 1868 1916];
ppos = [186 244 297 311 334 410 417 441 449 760 769 806 813 1001 1027 1030 1048 1105 1146 1176 1203 1350 1380 1470 1604 1606 1610 1624 1731 1898 1954 2016]; % �tima distribui��o de frequ�ncia (subcarriers) para 32 pilotos.

fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:M-1)/M);               % Normalized pilot signal, i.e., unit power.
%P = 1/sqrt(2)*(randn(1,M) + 1i*randn(1,M));

Pdiag = diag(P);

P = P.';

%% ********* Fourier Basis *********
F = fft(eye(N));
Fl = F(ppos,:);

%% ********* Channel ************
h0 = 0.5 + 1i*0.8;
h1 = 0.7 - 1i*0.3;
h2 = 0.1 - 1i*0.4;

energy = (abs(h0).^2 + abs(h1).^2 + abs(h2).^2) / 3;
normalization_factor = 1/sqrt(energy);

h0 = h0*normalization_factor;
h1 = h1*normalization_factor;
h2 = h2*normalization_factor;

channel_gains = [h0, h1, h2];

h = complex(zeros(N,1),zeros(N,1));
pos = [1 24 31];
L = pos(length(pos));
for idx=1:1:length(pos)
    h(pos(idx)) = channel_gains(idx);
end

channel_energy = sum(abs(h).^2)/length(pos);

K = length(pos);                                % Number of non-zero entries

% Comments: The channel "h" has a structured model. You can insert more realistics channel models.
H = Fl*h; % H is the FFT of the channel H computed at the pilot's frequencies.

% Hl = fft(h,N);
% Hl = Hl(ppos,:);
% H_error = sum(abs(H-Hl).^2)/length(H);

%% ************* Transmission ****************

% Transmitted signal in frequency domain.
s = Pdiag*H;

% noise
rng(839);

SNR = -4; % SNR given in dB.
linearSNR = 10^(-SNR/20);
noise = linearSNR*((randn(size(s)) + 1i*randn(size(s))) / sqrt(2));

% received signal computed at pilot's frequencies.
%y = s + noise;

y = s;

%% ************ L1 Convex Estimation (CE) Form 1 ***************
% Sparse Channel estimation
pathGains = [h0 h1 h2];
lambda1 = M*10^(-SNR/10)/sum(abs(pathGains));

A = Pdiag*Fl(:,1:pos(length(pos)));
cvx_begin quiet
variable x(L) complex
% sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
%minimize( quad_form(A*x-y,eye(M))+lambda1*norm(x,1) )
minimize( quad_form(A*x-y,eye(M))+lambda1*norm(x,1) )
cvx_end
% building channel at all locations (simply from the dictionary).
%H_Sparse = Fl*x;

g_hat_ce1 = x;

error_ce1 = norm(h(1:length(g_hat_ce1))-g_hat_ce1)^2/norm(h(1:length(g_hat_ce1)))^2;

fprintf(1,'error_ce1: %d\n',error_ce1);

%% ************ Compressed Sensing Estimation (CSE) ************

A = Pdiag*Fl(:,1:pos(length(pos)));

g_hat_omp = OMP_orig(A,y,K);

% Comments: The performance of OMP depends on the choice of F and P. In
% this example, I assume that the channel has Fourier structure. In more
% realistic channel problem, Fourier basis presents a power leakage which
% increases the number of non-zero entries of g. Therefore, the choice of
% the F impacts directly on the algorithm performance.

error_omp = norm(h(1:length(g_hat_omp))-g_hat_omp)^2/norm(h(1:length(g_hat_omp)))^2;

fprintf(1,'error_omp: %d\n',error_omp);

%% ************ L1 Convex Estimation (CE) Form 2 ***************
epsilon = 0.000000000002;

% Sparse Channel estimation
A = Pdiag*Fl(:,1:NCP);
cvx_begin quiet
variable x(NCP) complex
% sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
minimize( norm(x,1) )
subject to
norm(A*x-y,2)<=epsilon
cvx_end

g_hat_ce2 = x;

% Error
error_ce2 = norm(h(1:length(g_hat_ce2))-g_hat_ce2)^2/norm(h(1:length(g_hat_ce2)))^2;

fprintf(1,'error_ce2: %d\n',error_ce2);

%% ************ BDPN ***************
% In the presence of noisy or imperfect data, however, it is undesirable to exactly fit the
% linear system. Instead, the constraint in (BP) is relaxed to obtain the basis pursuit
% denoise (BPDN) problem
opts = spgSetParms('verbosity',0);

% Compressed Sensing
A = Pdiag*Fl(:,1:NCP);

% Reference: https://www.math.ucdavis.edu/~mpf/spgl1/index.html
% If there was noise inherent in the signal that we measured, it would be better to instead solve the basis pursuit denoise (BPDN) problem
% minimize ||x||_1  subject to  ||Ax - b||_2 <= sigma,

sigma = 0.000000000002; %linearSNR; %      % Desired ||Ax - b||_2
g_hat_bdpn = spg_bpdn(A, y, sigma, opts);

% Error
error_bdpn = norm(h(1:length(g_hat_bdpn))-g_hat_bdpn)^2/norm(h(1:length(g_hat_bdpn)))^2;

fprintf(1,'error_bdpn: %d\n',error_bdpn);

%% ************ CoSamp 2 ***************
A = Pdiag*Fl(:,1:pos(length(pos)));

%valores otimos
tol = 0.01;
maxiterations = 1000;
g_hat_cosamp2 = cosamp_mod(A,y,K,tol,maxiterations);

% Error
error_cosamp2 = norm(h(1:length(g_hat_cosamp2))-g_hat_cosamp2)^2/norm(h(1:length(g_hat_cosamp2)))^2;

fprintf(1,'error_cosamp2: %d\n',error_cosamp2);

%% ************ OMP 2 ***************
A = Pdiag*Fl(:,1:pos(length(pos)));

%valores otimos
tol = 0.01;
maxiterations = 1000;
g_hat_omp2 = cosaomp(A,y,K,tol,maxiterations);

% Error
error_omp2 = norm(h(1:length(g_hat_omp2))-g_hat_omp2)^2/norm(h(1:length(g_hat_omp2)))^2;

fprintf(1,'error_omp2: %d\n',error_omp2);


%% ************ CoSamp 3 ***************
A = Pdiag*Fl(:,1:pos(length(pos)));

%valores otimos
Its = 1000;
[g_hat_omp3,xcosamp] = cosamp_rice(y, A, K, Its);

% Error
error_omp3 = norm(h(1:length(g_hat_omp3))-g_hat_omp3)^2/norm(h(1:length(g_hat_omp3)))^2;

fprintf(1,'error_cosamp3: %d\n',error_omp3);

%% ************ Orthogonal Matching Pursuit 3 ***************
A = Pdiag*Fl(:,1:pos(length(pos)));


%yfit = wmpalg('BMP',y,A,'itermax',100);
%[ws,r] = temporalMP(y,A);
[g_hat_omp4, numIts] = romp(K, A, y);

% Error
error_mp1 = norm(h(1:length(g_hat_omp4))-g_hat_omp4)^2/norm(h(1:length(g_hat_omp4)))^2;

fprintf(1,'error_omp4: %d\n',error_mp1);


%% ************ Matching Pursuit 1 ***************
A = Pdiag*Fl(:,1:pos(length(pos)));

[g_hat_mp, tap_pos, res] = mpv2(A,y,K);
%[X_mp, E_mp] = mpv2(A,y);
%[s, err_mse, iter_time]=greed_mp(y,A,31);
%YFIT = wmpalg('BMP',y,A);

% Error
error_mp1 = norm(h(1:length(g_hat_mp))-g_hat_mp)^2/norm(h(1:length(g_hat_mp)))^2;

fprintf(1,'error_mp1: %d\n',error_mp1);


