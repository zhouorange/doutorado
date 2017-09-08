clear all;close all;clc

rng(1)

SNR = -10:2:18;                     % Signal-to-noise ratio in dB.

M = 30;                             % Number of antennas.
K = 10;                             % Number of single-antenna users.
P = 20;                             % Channel Length (Finite Impulse Response - FIR).
L = 7;                              % Number of cells.

N = getPilotLength(K,P);            % Pilot length is set according to K and P.
q = 10.^(SNR./10);                  % Uplink pilot power.

cellRadius = 1000;
cellHole = 100;
sshadow = 8;
gamma = 3.8;

NUM_ITER_BETA = 10000;

% Pilot signal generation
S = [];
for k_idx=1:1:K
    s_aux = generatePilotMatrixv1(N, P, k_idx); % No futuro deve-se usar a versão 2, pois ela normaliza a energia dos pilotos.
    S = [S, s_aux];
end

theoretical_ls_error_vec = zeros(1,length(q));
theoretical_mmse_error_vec = zeros(1,length(q));
prop_error_vec = zeros(1,length(q));
theoretical_prop_error_vec = zeros(1,length(q));
theoretical_prop_error_hat1_vec = zeros(1,length(q));
theoretical_prop_error_hat2_vec = zeros(1,length(q));
for q_idx=1:1:length(q)
    
    theoretical_ls_error = 0;
    theoretical_mmse_error = 0;
    theoretical_prop_error = 0;
    prop_error = complex(zeros(P,P),zeros(P,P));
    theoretical_prop_error_hat1 = 0;
    theoretical_prop_error_hat2 = 0;
    for beta_iter=1:1:NUM_ITER_BETA
        
        % Generate Beta (large scale coefficients).
        beta = GenerateRandomBetav4(L,K,cellRadius,cellHole,sshadow,gamma);
        
        % Received pilot symbols from all users in all the L cells at
        % BS number 1 that is the cell at the middle, i.e., with L-1
        % cells around it.
        Y1 = complex(zeros(M,N),zeros(M,N));
        for l_idx=1:1:L
            
            G = [];
            for k_idx=1:1:K
                h = (1/sqrt(2)).*complex(randn(M,P),randn(M,P));
                g = sqrt(beta(l_idx,k_idx)).*h;
                G = [G, g];
            end
            
            if(l_idx==1) % In our simulations we select only cell 1.
                G11_all_k = G;
            end
            
            Y1 = Y1 + sqrt(q(q_idx))*G*(S');
            
        end
        
        % Noise generation
        W = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Add noise to received signal.
        Y1 = Y1 + W;
        
        % *************************** Channel estimation ***************************
        
        for k_idx=1:1:K
            
            var1 = 0.1;
            var2 = 0.01;

            [beta_sum_ik, beta11k, beta11k_hat1, beta11k_hat2, beta_sum_ik_hat1, beta_sum_ik_hat2] = getBetaSumv2(L, beta, k_idx, var1, var2);
            
            epsilon1k = (beta_sum_ik + 1/(q(q_idx)*N));
            
            epsilon1k_hat1 = (beta_sum_ik_hat1 + 1/(q(q_idx)*N));
            
            epsilon1k_hat2 = (beta_sum_ik_hat2 + 1/(q(q_idx)*N));
            
            % Select the channel
            G11k = G11_all_k(:,(P*(k_idx - 1) + 1) : (P*k_idx));
            
            % Select only the pilots of the k-th user.
            Sk = S(:,(P*(k_idx - 1) + 1) : (P*k_idx));
            
            % 1) Least Squares (LS) Solution (Estimation).
            Z11_ls = (1/(sqrt(q(q_idx))*N))*Y1*Sk;
            
            % Auxiliar variable.
            aux = ((Z11_ls')*(Z11_ls));
            
            % 3) Proposed estimation.
            Z11_prop = (M*beta11k*Z11_ls)./(trace(aux)./P);
            
            prop_error = prop_error + ((Z11_prop'-G11k')*(Z11_prop-G11k));
            
            % Accumulate LS error for average.
            theoretical_ls_error = theoretical_ls_error + (epsilon1k - beta11k);
            % Accumulate MMSE error for average.
            theoretical_mmse_error = theoretical_mmse_error + ((beta11k/epsilon1k) * (epsilon1k - beta11k));
            % Accumulate Prop. error for average.
            theoretical_prop_error = theoretical_prop_error + (beta11k*((((2-M*P)*beta11k)/((P*M-1)*epsilon1k)) + 1));
  
            % Accumulate Prop. error for average.
            theoretical_prop_error_hat1 = theoretical_prop_error_hat1 + (beta11k_hat1*((((2-M*P)*beta11k_hat1)/((P*M-1)*epsilon1k_hat1)) + 1));       
            % Accumulate Prop. error for average.
            theoretical_prop_error_hat2 = theoretical_prop_error_hat2 + (beta11k_hat2*((((2-M*P)*beta11k_hat2)/((P*M-1)*epsilon1k_hat2)) + 1));
            
        end
        
    end
    
    fprintf('SNR: %d\n',SNR(q_idx));
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error_vec(q_idx) = theoretical_ls_error/(NUM_ITER_BETA*K);
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error_vec(q_idx) = theoretical_mmse_error/(NUM_ITER_BETA*K);
    
    % Proposed Estimator Simulation Error. (Monte Carlo)
    prop_error = prop_error./(M * NUM_ITER_BETA * K);
    prop_error_vec(q_idx) = sum(diag(prop_error))/P;
    
    % Proposed Estimator Theoretical error.
    theoretical_prop_error_vec(q_idx) = theoretical_prop_error./(NUM_ITER_BETA*K);
    
    % Proposed Estimator Simulation Error - sigma^2 = 0.1. (Monte Carlo)
    theoretical_prop_error_hat1_vec(q_idx) = theoretical_prop_error_hat1./(NUM_ITER_BETA*K);
    
    % Proposed Estimator Simulation Error - sigma^2 = 0.01. (Monte Carlo)
    theoretical_prop_error_hat2_vec(q_idx) = theoretical_prop_error_hat2./(NUM_ITER_BETA*K);
    
end

fdee_figure = figure;
semilogy(SNR,real(theoretical_mmse_error_vec),'--r');
hold on;
semilogy(SNR,real(theoretical_ls_error_vec),'--b','MarkerSize',7);
semilogy(SNR,real(theoretical_prop_error_vec),'--c','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'c*','MarkerSize',7);
semilogy(SNR,real(theoretical_prop_error_hat1_vec),'ko','MarkerSize',7);
semilogy(SNR,real(theoretical_prop_error_hat2_vec),'k*','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d',M,K,P,L,N);
title(titleStr);
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','LS (ana)','Prop. (ana)','Prop. (sim)','Prop. (sim) \sigma^{2}=0.1','Prop. (sim) \sigma^{2}=0.01');

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v24_%s.fig',M,K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v24_%s.mat',M,K,P,L,N,timeStamp);
save(fileName);
