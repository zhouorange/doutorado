clear all;close all;clc

rng(1)

SNR = 10;           % Signal-to-noise ratio in dB.

M = 30;             % Number of antennas.
K = 10;             % Number of single-antenna users.
P = 20;             % Channel Length (Finite Impulse Response - FIR).
L = 7;              % Number of cells.

N = getPilotLength(K,P);  % Pilot length is set according to K and P.
q = 10.^(SNR./10);            % Uplink pilot power.

isFixedValue = true;
if(isFixedValue)
    a = 0.05;
    beta = a*ones(L,K);
    for ll=1:1:L
        for kk=1:1:K
            if(kk==ll)
                beta(ll,kk) = 1;
            end
        end
    end
else
    beta = abs(randn(L,K));
end
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);

q_idx=1;

epsilon11 = (beta_sum + 1/(q(q_idx)*N));

numIter = 10000;
trace_res_sqr = 0;
var_b = zeros(1,numIter);
est_g_prop_vs_g = zeros(P,P);
for kk=1:1:numIter
    
    % Pilot signal generation
    S = [];
    for k_idx=1:1:K
        s_aux = generatePilotMatrixv1(N, P, k_idx);
        %s_aux = generatePilotMatrixFFT(N, k_idx);
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
        
        Y1 = Y1 + sqrt(q(q_idx))*C*(S');
        
    end
    
    % Add noise to received signal.
    Y1 = Y1 + W;
    
    % 1) Least Squares (LS) Solution (Estimation).
    Z11_ls = (1/(sqrt(q(q_idx))*N))*Y1*S1;
        
    sqr_matrix = Z11_ls'*Z11_ls;
    
    Z11_prop = (M*P*beta111*Z11_ls)./(trace(sqr_matrix));
    
    est_g_prop_vs_g = est_g_prop_vs_g + (Z11_prop'*C111);
    
end

estimated_value = trace(est_g_prop_vs_g/numIter)/P;

correct_value = M*(beta111.^2)/epsilon11;

fprintf(1,'estimated_value: %d\n',estimated_value);
fprintf(1,'correct_value: %d\n',correct_value);
