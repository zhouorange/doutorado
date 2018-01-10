clear all;close all;clc

rng(1)

M = 30;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.

a = 0.05;           % Constant beta value.
beta111 = 1;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

% Simulation loop starts here.
NUM_ITER = 100000;

for numIter = 1:1:NUM_ITER
    
    beta_sum = 0;
    sum_G = zeros(M,N);
    sum_aux_G = zeros(M,N);
    Gil = zeros(M,K,L);
    % Iterate over all cells (L) in the assumed system.
    % Here we consider the target cell, i, is 1, i.e., i = 1.
    for l=1:1:L
        
        % Generate channel matrix G_{il}.
        beta = a;
        if(l == 1)
            beta = beta111;
        end
        betaMatrix = sqrt(beta)*eye(K);
        Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
        
        % Summation of all channels.
        sum_G = sum_G + (Gil(:,:,l)*(S'));
        
        sum_aux_G = sum_aux_G + Gil(:,:,l);
        
    end
    
    sum_aux_G = sum_aux_G*(S');
    
    error = sum_G - sum_aux_G;
    
end


