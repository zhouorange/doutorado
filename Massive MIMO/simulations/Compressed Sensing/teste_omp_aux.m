clear all;
clc;
%%
% MIMO Channel 
N = 100;              % number of antennas
F = fft(eye(N));        % Fourier Basis

h0 = 0.5 + 1i*0.8;
h1 = 0.7 - 1i*0.3;
h2 = 0.1 - 1i*0.4;

energy = (abs(h0).^2 + abs(h1).^2 + abs(h2).^2) / 3;
normalization_factor = 1/sqrt(energy);

h0 = h0*normalization_factor;
h1 = h1*normalization_factor;
h2 = h2*normalization_factor;

g = zeros(N,1);
pos = [1 24 31];
g(pos) = [h0 h1 h2]; % Sparse representation

h= F*g;

K=3;  % Number of non-zero entries

% Comments: The channel "h" has a structured model. You can
% insert more realistics channel models.
     


% Pilot sequence 
M=12;                                            % Length of the training sequence
P=1/sqrt(2)*(randn(M,N) + 1i*randn(M,N));
% Comments : The matix P addresses pilot sequences to each antenna. 

s=P*h;                           % Transmit signal

% noise
SNR = 10; % SNR given in dB.
linearSNR = 10^(-SNR/20);
noise = linearSNR*((randn(size(s)) + 1i*randn(size(s))) / sqrt(2));


% received signal
x = s + noise;


% Estimation of g
[g_est]=omp(P*F,x,K);

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

error = norm(g-g_est)^2/norm(g)^2

%email: araujo@gtel.ufc.br / daniel.c.araujo@gmail.com

