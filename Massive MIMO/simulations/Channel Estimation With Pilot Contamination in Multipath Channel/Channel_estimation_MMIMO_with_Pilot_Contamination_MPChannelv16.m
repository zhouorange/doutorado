clear all;close all;clc

rng(1)

SNR = 10;                           % Signal-to-noise ratio in dB.

M = 600;                            % Number of antennas.
K = 10;                             % Number of single-antenna users.
P = 20;                             % Channel Length (Finite Impulse Response - FIR).
L = 7;                              % Number of cells.

N = getPilotLength(K,P);            % Pilot length is set according to K and P.
q = 10.^(SNR./10);                  % Uplink pilot power.

a = logspace(-3, 1, 20);

NUM_ITER = 10;

isFixedValue = true;

theoretical_ls_error = zeros(1,length(a));
theoretical_mmse_error = zeros(1,length(a));
prop_4_ana_error_vec = zeros(1,length(a));
for a_idx=1:1:length(a)
    
    % Calculate Beta.
    [beta_sum, beta111, beta] = generateBetav0(a(a_idx), L, K, isFixedValue);
    
    for q_idx=1:1:length(q)
        
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
                
        fprintf('a: %d\n',a(a_idx));
        
        % Least Squares (LS) Theoretical error.
        theoretical_ls_error(a_idx) = epsilon11 - beta111;
        
        % Minimum Mean Squared Error (MMSE) Theoretical error.
        theoretical_mmse_error(a_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
        
        % Proposed Estimator Theoretical error.
        prop_4_ana_error_vec(a_idx) = beta111*( ( ((2-M*P)*beta111)/((P*M-1)*epsilon11)  ) + 1);
        
    end
end

fdee_figure = figure;
loglog(a,theoretical_mmse_error,'--r');
hold on;
loglog(a,theoretical_ls_error,'--b','MarkerSize',7);
loglog(a,real(prop_4_ana_error_vec),'c*','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d',M,K,P,L,N);
title(titleStr);
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)');

% % Get timestamp for saving files.
% timeStamp = datestr(now,30);
% fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v16_%s.fig',M,K,P,L,N,timeStamp);
% savefig(fdee_figure,fileName);
% 
% % Save workspace to .mat file.
% clear fdee_figure
% fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v16_%s.mat',M,K,P,L,N,timeStamp);
% save(fileName);
