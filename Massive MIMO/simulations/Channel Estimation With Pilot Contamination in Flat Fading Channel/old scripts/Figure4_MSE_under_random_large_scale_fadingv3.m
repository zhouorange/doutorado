clear all;close all;clc

rng(1)

SNR = -10:4:18;     % Signal-to-noise ratio in dB.
linear_SNR = 10.^(SNR./10);

M = 30;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.

cellRadius = 1000;
cellHole = 100;
sshadow = 8;
gamma = 3.8;

NUM_ITER = 1000;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

theoretical_ls_error = zeros(1,length(SNR));
theoretical_mmse_error = zeros(1,length(SNR));
theoretical_proposed_error = zeros(1,length(SNR));
theoretical_proposed_approx_error = zeros(1,length(SNR));
ls_error_vec = zeros(1,length(SNR));
mmse_error_vec = zeros(1,length(SNR));
prop_error_vec = zeros(1,length(SNR));
for snr_idx=1:1:length(SNR)
    
    for numIter = 1:1:NUM_ITER
        
        % Generate Beta (large scale coefficients).
        beta = GenerateRandomBetav4(L,K,cellRadius,cellHole,sshadow,gamma);
        
        sum_of_all_betas = getSumOfAllBetas(beta,L,K);
        
        % Uplink pilot power.
        q = (K*linear_SNR(snr_idx))./(N*(sum_of_all_betas)); 
        
        sum_G = zeros(M,N);
        Gil = zeros(M,K,L);
        for l=1:1:L
            
            % Generate channels.
            betaMatrix = diag(sqrt(beta(l,:)));
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
            
            % Summation of all channels.
            sum_G = sum_G + Gil(:,:,l)*(S');
            
        end
        
        % Apply squared pilot power.
        sum_G = sqrt(q)*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % received pilot symbols at BS.
        Y1 = sum_G + W1;
        
        %-------------- Channel Estimation --------------
        
        % Summation of betas.
        [beta_sum, beta111] = getBetaSum(beta, 1, L);
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q*N));
        
        % ******* LS ********
        Z11_ls = (1/(sqrt(q)*N))*Y1*S(:,1);
        ls_error_vec(snr_idx) = ls_error_vec(snr_idx) + ((Z11_ls'-Gil(:,1,1)')*(Z11_ls-Gil(:,1,1)));
        
        % Least Squares (LS) Theoretical error.
        theoretical_ls_error(snr_idx) = theoretical_ls_error(snr_idx) + (epsilon11 - beta111);
        
        % ******* MMSE ********
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        mmse_error_vec(snr_idx) = mmse_error_vec(snr_idx) + ((Z11_mmse'-Gil(:,1,1)')*(Z11_mmse-Gil(:,1,1)));
        
        % Minimum Mean Squared Error (MMSE) Theoretical error.
        theoretical_mmse_error(snr_idx) = theoretical_mmse_error(snr_idx) + ((beta111/epsilon11) * (epsilon11 - beta111));
        
        % ******* Prop. ********
        Z11_prop = M*beta111*Z11_ls/((Z11_ls')*(Z11_ls));
        prop_error_vec(snr_idx) = prop_error_vec(snr_idx) + ((Z11_prop'-Gil(:,1,1)')*(Z11_prop-Gil(:,1,1)));
        
        % Proposed estimator Analytical error.
        theta_ik = calculateTheta_ikv2(beta111, epsilon11, M);
        theoretical_proposed_error(snr_idx) = theoretical_proposed_error(snr_idx) + (((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik);        
        
        % Proposed estimator Approximated error.
        theoretical_proposed_approx_error(snr_idx) = theoretical_proposed_approx_error(snr_idx) + ((beta111*((((2-M)*beta111)/((M-1)*epsilon11)) + 1)));
    end
    
    fprintf('SNR: %d\n',SNR(snr_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error_vec(snr_idx) = ls_error_vec(snr_idx)./(M * NUM_ITER);
    % Least Squares (LS) Analytical Error.
    theoretical_ls_error(snr_idx) = theoretical_ls_error(snr_idx)./NUM_ITER;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error_vec(snr_idx) = mmse_error_vec(snr_idx)./(M * NUM_ITER);
    % Minimum Mean Squared Error (MMSE) Analytical Error.
    theoretical_mmse_error(snr_idx) = theoretical_mmse_error(snr_idx)./NUM_ITER;
    
    % Proposed estimator Simulation Error (Monte Carlo)
    prop_error_vec(snr_idx) = prop_error_vec(snr_idx)./(M * NUM_ITER);
    
    % Proposed estimator Analytical error.
    theoretical_proposed_error(snr_idx) = theoretical_proposed_error(snr_idx)./NUM_ITER;
    
    % Proposed estimator Approximated Error
    theoretical_proposed_approx_error(snr_idx) = theoretical_proposed_approx_error(snr_idx)./NUM_ITER;
end

fontSize = 10;

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
ylabel('avg. MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (sim)', 'Prop. (approx.)');
%axis([SNR(1) SNR(length(SNR)) 0.18 5.1])
strText = sprintf('M = %d',M);
x1 = SNR(length(SNR))-12;
y1 = 1;
text(x1,y1,strText,'FontSize', fontSize)

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure4_MSE_under_random_large_scale_fading_M%d_K%d_L%d_N%d_v0_%s.fig',M,K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure4_MSE_under_random_large_scale_fading_M%d_K%d_L%d_N%d_v0_%s.mat',M,K,L,N,timeStamp);
save(fileName);
