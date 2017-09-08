clear all;close all;clc

rng(1)

SNR = -20:4:22;     % Signal-to-noise ratio in dB.

M = 30;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.
q = 10.^(SNR./10);  % Uplink pilot power.

a = 0.05;           % Constant beta value.
beta111 = 1;

NUM_ITER = 10000;

% Generate pilot signals.
S = generatePilotMatrixFFT(N);

theoretical_ls_error = zeros(1,length(q));
theoretical_mmse_error = zeros(1,length(q));
theoretical_proposed_error = zeros(1,length(q));
theoretical_proposed_approx_error = zeros(1,length(q));
ls_error_vec = zeros(1,length(q));
mmse_error_vec = zeros(1,length(q));
prop_error_vec = zeros(1,length(q));
for q_idx=1:1:length(q)
    
    ls_error = 0;
    mmse_error = 0;
    prop_error = 0;
    for numIter = 1:1:NUM_ITER
        
        beta_sum = 0;
        sum_G = zeros(M,K);
        Gil = zeros(M,K,L);
        for l=1:1:L
            
            % Generate channels.
            beta = a;
            if(l == 1)
                beta = 1;
            end
            betaMatrix = sqrt(beta)*eye(K);
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
            
            % Summation of all channels.
            sum_G = sum_G + Gil(:,:,l)*(S');
            
            % Summation of betas.
            beta_sum = beta_sum + beta;
            
        end
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,K),randn(M,K));
        
        % received pilot symbols at BS.
        Y1 = sum_G + W1;
        
        % ******* LS ********
        Z11_ls = (1/(sqrt(q(q_idx))*N))*Y1*S(:,1);       
        ls_error = ls_error + ((Z11_ls'-Gil(:,1,1)')*(Z11_ls-Gil(:,1,1)));
        
        % ******* MMSE ********
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        mmse_error = mmse_error + ((Z11_mmse'-Gil(:,1,1)')*(Z11_mmse-Gil(:,1,1)));
        
        % ******* Prop. ********
        Z11_prop = M*beta111*Z11_ls/((Z11_ls')*(Z11_ls));
        prop_error = prop_error + ((Z11_prop'-Gil(:,1,1)')*(Z11_prop-Gil(:,1,1)));
    end
    
    fprintf('SNR: %d\n',SNR(q_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M * NUM_ITER);
    ls_error_vec(q_idx) = ls_error;
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error(q_idx) = epsilon11 - beta111;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M * NUM_ITER);
    mmse_error_vec(q_idx) = mmse_error;
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error(q_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
    
    % Proposed estimator Simulation Error (Monte Carlo)
    prop_error = prop_error./(M * NUM_ITER);
    prop_error_vec(q_idx) = prop_error;   
    
    % Proposed estimator Theoretical error.
    theta_ik = calculateTheta_ikv1(beta111, epsilon11, M);
    %theta_ik = 0.354626695759462; % para valores de M > 85
    theoretical_proposed_error(q_idx) = ((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
    
    % Proposed estimator Approximated error.
    theoretical_proposed_approx_error(q_idx) = (beta111*((((2-M)*beta111)/((M-1)*epsilon11)) + 1));
end

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'--r');
hold on;
semilogy(SNR,real(mmse_error_vec),'r*','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'--b','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b*','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_error),'--k','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'k*','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_approx_error),'ko','MarkerSize',7);
hold off
grid on;
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (sim)', 'Prop. (approx.)');

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure1_MSE_versus_TX_SNR_M%d_K%d_L%d_N%d_v0_%s.fig',M,K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure1_MSE_versus_TX_SNR_M%d_K%d_L%d_N%d_v0_%s.mat',M,K,L,N,timeStamp);
save(fileName);
