clear all;close all;clc

rng(1)

SNR = 100;          % Signal-to-noise ratio in dB.

M = 30;            % Number of antennas.
K = 2;              % Number of single-antenna users.
P = 4;              % Channel Length (Finite Impulse Response - FIR).
L = 2;              % Number of cells.

N = getPilotLength(K,P);  % Pilot length is set according to K and P.
q = 100000000;            % Uplink pilot power.

NUM_ITER = 10000;

beta = abs(randn(L,K));
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);
epsilon11 = (beta_sum + 1/(q*N));

ls_error = complex(zeros(P,P),zeros(P,P));
mmse_error = complex(zeros(P,P),zeros(P,P));
prop_error = complex(zeros(P,P),zeros(P,P)); 
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
    Z11_mmse = (beta111/epsilon11)*Z11_ls;
    
    mmse_error = mmse_error + ((Z11_mmse'-C111')*(Z11_mmse-C111));
    
    % 3) Proposed Estimation.
    
    aux = (Z11_ls')*(Z11_ls);
    
    Z11_prop = M*beta111*(Z11_ls * ((aux)^(-1)));
    
    prop_error = prop_error + ((Z11_prop'-C111')*(Z11_prop-C111));
    
end

% Least Squares (LS) Simulation Error. (Monte Carlo)
ls_error = ls_error./(M * NUM_ITER);

% Least Squares (LS) Theoretical error.
theoretical_ls_error = (beta_sum + 1/(q*N)) - beta111;

% Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
mmse_error = mmse_error./(M * NUM_ITER);

% Minimum Mean Squared Error (MMSE) Theoretical error.
theoretical_mmse_error = (beta111/epsilon11) * (epsilon11 - beta111);

% Proposed Estimator Simulation Error. (Monte Carlo)
prop_error = prop_error./(M * NUM_ITER);

theta_ik = calculateTheta_ik(beta111, epsilon11, M);
theoretical_proposed_error = ((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;


