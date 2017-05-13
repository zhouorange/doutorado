clear all;
clc;
%%
% MIMO Channel 
N=100;              % number of antennas
F=fft(eye(N));        % Fourier Basis

%%%%% DEBUG %%%%%%%%%%
g = zeros(N,1);
pos = [1 24 31];
g(pos) = [(1.2+1i*3.4) (5-1i*2) (2.1+1i*1.1)];
%%%%% DEBUG %%%%%%%%%%

h= F*g;

%%%%% DEBUG %%%%%%%%%%
K = length(pos); % sum(find(g>=0.01));  % Number of non-zero entries
%%%%% DEBUG %%%%%%%%%%

% Comments: The channel "h" has a structured model. You can
% insert more realistics channel models.
     


% Pilot sequence 
M=ceil(1.001*K);                          % Length of the training sequence
P=1/sqrt(2)*(randn(M,N) + 1i*randn(M,N));
% Comments : The matix P addresses pilot sequences to each antenna. 

s=P*h;                           % Transmit signal

% noise
sigma=0.1;
noise= sqrt(sigma/2)*(randn(size(s))+1i*randn(size(s)));


% received signal
x = s + noise;

%x = s;


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

