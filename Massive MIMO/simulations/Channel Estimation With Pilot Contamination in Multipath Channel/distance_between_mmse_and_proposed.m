clear all;close all;clc

rng(1)

SNR = 10;             % Signal-to-noise ratio in dB.

M = 30;%30; %85;            % Number of antennas.
K = 10;                     % Number of single-antenna users.
P = 20;                     % Channel Length (Finite Impulse Response - FIR).
L = 7;                      % Number of cells.

N = getPilotLength(K,P);    % Pilot length is set according to K and P.
q = 10.^(SNR./10);          % Uplink pilot power.

NUM_ITER = 10000;

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

ls_error_vec = zeros(1,length(q));
theoretical_ls_error = zeros(1,length(q));
mmse_error_vec = zeros(1,length(q));
theoretical_mmse_error = zeros(1,length(q));
prop_error_vec = zeros(1,length(q));
theoretical_distance = zeros(1,length(q));
for q_idx=1:1:length(q)
    
    epsilon11 = (beta_sum + 1/(q(q_idx)*N));
    
    distance = zeros(P,P);
    for iter=1:1:NUM_ITER
        
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
        
        % 2) Minimum Mean Squared Error (MMSE) Solution (Estimation).
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        
        % 3) Proposed Estimator.
        aux = ((Z11_ls')*(Z11_ls));
        Z11_prop = (M*P*beta111*Z11_ls)./(trace(aux));
        
        distance = distance + ((Z11_prop'-Z11_mmse')*(Z11_prop-Z11_mmse));
        
    end
    
    fprintf('SNR: %d\n',SNR(q_idx));
    
    %----------------------- distance between Proposed and MMSE estimators -----------------------
    % Proposed Estimator Simulation Error. (Monte Carlo)
    distance = distance./(M*NUM_ITER);
    prop_error_vec(q_idx) = sum(diag(distance))/P;
    
    % Proposed Estimator Theoretical error.
    theoretical_distance(q_idx) = (1/(P*M-1))*((beta111.^2)./epsilon11);
    
    fprintf(1,'Estimated distance: %d\n',prop_error_vec(q_idx));
    fprintf(1,'Theoretical distance: %d\n',theoretical_distance(q_idx));
end

