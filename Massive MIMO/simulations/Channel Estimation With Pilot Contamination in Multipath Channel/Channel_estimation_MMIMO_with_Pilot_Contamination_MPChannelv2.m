clear all;close all;clc

SNR = 100;          % Signal-to-noise ratio in dB.

M = 500;            % Number of antennas.
K = 2;              % Number of single-antenna users.
P = 4;              % Channel Length (Finite Impulse Response - FIR).
L = 2;              % Number of cells.

N = 13;             % Pilot length.
q = 100000000;            % Uplink pilot power.

NUM_ITER = 100000;

beta = abs(randn(L,K));

ls_error = complex(zeros(P,P),zeros(P,P));
theoretical_ls_error = 0;
for iter=1:1:NUM_ITER
    
    % Pilot signal generation
    S = [];
    for k_idx=1:1:K
        s_aux = generatePilotMatrixv1(N, P, k_idx);
        if(k_idx==1)
            S1 = s_aux;
        end
        S = [S, s_aux];
    end
    
    % Noise generation
    W = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
    
    % Received pilot symbols at BS i.
    Y1 = complex(zeros(M,N),zeros(M,N));
    for l_idx=1:1:L
        
        C = [];
        for k_idx=1:1:K
            g = (1/sqrt(2))*complex(randn(M,P),randn(M,P));
            c = sqrt(beta(l_idx,k_idx))*g;
            C = [C, c];
            
            if(l_idx==1 && k_idx==1)
                C111 = c;
            end
        end
        
        Y1 = Y1 + sqrt(q)*C*(S');
        
    end
    
    % Add noise to received signal.
    Y1 = Y1 + W;
    
    % 1) Least Squares (LS) Solution (Estimation).
    Z11_ls = (1/(sqrt(q)*N))*Y1*S1;

    ls_error = ls_error + ((Z11_ls'-C111')*(Z11_ls-C111));
    
    % 2) Minimum Mean Squared Error (MMSE) Solution (Estimation).
    Z11_mmse = 1;
    
end

% Simulation Error. (Monte Carlo)
ls_error = ls_error./(M * NUM_ITER);

% Least Squares (LS) Theoretical error.
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);
theoretical_ls_error = theoretical_ls_error + (beta_sum + 1/(q*N)) - beta111;
