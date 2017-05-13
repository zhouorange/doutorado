clear all;close all;clc;

N = 512;                % number of OFDM subcarriers.
F = fft(eye(N));        % Fourier Basis

% Comments: The vector P addresses pilot sequences.
% Pilot sequence 
M = 32;                                         % Number of pilots. Length of the training sequence
pos = 1:N/M:N;                                  % Position of the pilot in the frequency domain.
%P = 1/sqrt(2)*(randn(1,M) + 1i*randn(1,M));

fo = 1000; % in Hz.
P = exp(1i*2*pi*fo*(0:1:M-1)/M)/sqrt(M);

Pdiag = diag(P);

Fl = F(pos,:);

%%
% MIMO Channel 
%g = exp(-0.7*(0:N-1)'); % Sparse representation.






K = length(pos); % sum(find(g>=0.01));  % Number of non-zero entries

% Comments: The channel "h" has a structured model. You can
% insert more realistics channel models.
H = Fl*g; % h is the FFT of the channel g.

s = Pdiag*H;                % Transmit signal in frequency domain.

% noise
rng(839);
sigma = 0.1;



%noise = sqrt(sigma/2)*(randn(size(s))+1i*randn(size(s)));

SNR = 10;
linearSNR = 10^(SNR/10);
noise = (randn(size(s)) + 1i*randn(size(s))) / sqrt(2) / ( linearSNR );


% received signal
x = s + noise;

%x = s;

% Estimation of g
A = Pdiag*Fl;
%[g_est] = omp(A,x,K);
[g_est] = OMP_orig(A,x,K);

% Comments: The performance of OMP depends on the choice of F and P. In
% this example, I assume that the channel has Fourier structure. In more
% realistic channel problem, Fourier basis presents a power leakage which
% increases the number of non-zero entries of g. Therefore, the choice of
% the F impacts directly on the algorithm performance.

% You can add an input 'opt' in omp(A,x,K,opt). 

%   'opt'  is a structure with more options, including:
%       .target     = Define a residual in case you do not have knowladge
%                     of the sparsity.
%       .mode   = {0,1};  
%       If mode=0, there is no orthogonalization procedure. 
%       If mode=1 it is performed an orthonormalization on the matrix A.

[g g_est];

error_cse = norm(g-g_est)^2/norm(g)^2
%fprintf(1,'error: %f\n',error);

%email: araujo@gtel.ufc.br / daniel.c.araujo@gmail.com

%% ************* Least Squares Estimation (LSE) ******************

A = Pdiag*Fl(:,1:pos(length(pos)));

invMat = (((A'*A)^(-1))*A');

g_hat_lse = invMat*x;

error_lse = norm(g(1:length(g_hat_lse))-g_hat_lse)^2/norm(g(1:length(g_hat_lse)))^2

%% ************* Mean Squared Error Estimation (MMSE-E) ******************
invMat = (A*A' + (1/linearSNR)*eye(M))^-1;
invMat = A'*invMat;

g_hat_mmsee = invMat*x;

error_mmsee = norm(g(1:length(g_hat_mmsee))-g_hat_mmsee)^2/norm(g(1:length(g_hat_mmsee)))^2