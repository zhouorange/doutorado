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

k_index = 1;

NUM_ITER = 100;

ls_error_vec = zeros(1,length(q));
theoretical_ls_error = zeros(1,length(q));
mmse_error_vec = zeros(1,length(q));
theoretical_mmse_error = zeros(1,length(q));
prop_error_vec = zeros(1,length(q));
theoretical_prop_error = zeros(1,length(q));

% Calculate Beta.
[distances, angles, beta, beta11k, beta_sum] = GenerateRandomBetav3(L,K,cellRadius,cellHole,sshadow,gamma,k_index);

for q_idx=1:1:length(q)
    
    epsilon1k = (beta_sum + 1/(q(q_idx)*N));
    
    ls_error = complex(zeros(P,P),zeros(P,P));
    mmse_error = complex(zeros(P,P),zeros(P,P));
    prop_error = complex(zeros(P,P),zeros(P,P));
    for iter=1:1:NUM_ITER
        
        % Pilot signal generation
        S = [];
        for k_idx=1:1:K
            s_aux = generatePilotMatrixv1(N, P, k_idx);
            if(k_idx==k_index)
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
                g = (1/sqrt(2)).*complex(randn(M,P),randn(M,P));
                c = sqrt(beta(l_idx,k_idx)).*g;
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
        
        ls_error = ls_error + ((Z11_ls'-C111')*(Z11_ls-C111));
        
        % 2) Minimum Mean Squared Error (MMSE) Solution (Estimation).
        Z11_mmse = (beta11k/epsilon1k)*Z11_ls;
        
        mmse_error = mmse_error + ((Z11_mmse'-C111')*(Z11_mmse-C111));
        
        % Auxiliar variable.
        aux = ((Z11_ls')*(Z11_ls));
        
        % 3) Proposed estimation.
        Z11_prop = (M*beta11k*Z11_ls)./(trace(aux)./P);
        
        prop_error = prop_error + ((Z11_prop'-C111')*(Z11_prop-C111));
    end
    
    fprintf('SNR: %d\n',SNR(q_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M * NUM_ITER);
    ls_error_vec(q_idx) = sum(diag(ls_error))/P;
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error(q_idx) = epsilon1k - beta11k;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M * NUM_ITER);
    mmse_error_vec(q_idx) = sum(diag(mmse_error))/P;
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error(q_idx) = (beta11k/epsilon1k) * (epsilon1k - beta11k);
    
    % Proposed Estimator Simulation Error. (Monte Carlo)
    prop_error = prop_error./(M * NUM_ITER);
    prop_error_vec(q_idx) = sum(diag(prop_error))/P;
    
    % Proposed Estimator Theoretical error.
    theoretical_prop_error(q_idx) = beta11k*((((2-M*P)*beta11k)/((P*M-1)*epsilon1k)) + 1);
    
end

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'--r');
hold on;
semilogy(SNR,real(mmse_error_vec),'r*','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'--b','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b*','MarkerSize',7);
semilogy(SNR,theoretical_prop_error,'--c','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'c*','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d',M,K,P,L,N);
title(titleStr);
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)','LS (sim)','Prop. (ana)','Prop. (sim)');

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v18_%s.fig',M,K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v18_%s.mat',M,K,P,L,N,timeStamp);
save(fileName);
