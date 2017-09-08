clear all;close all;clc

rng(1)

SNR = 10;           % Signal-to-noise ratio in dB.

M = 10:10:200;      % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.
q = 10.^(SNR./10);  % Uplink pilot power.

a = 0.05;           % Interferance leval value.
beta111 = 1;

NUM_ITER = 100000;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

theoretical_ls_error = zeros(1,length(M));
theoretical_mmse_error = zeros(1,length(M));
theoretical_proposed_error = zeros(1,length(M));
theoretical_proposed_approx_error = zeros(1,length(M));
ls_error_vec = zeros(1,length(M));
mmse_error_vec = zeros(1,length(M));
prop_error_vec = zeros(1,length(M));
for m_idx=1:1:length(M)    
    ls_error = 0;
    mmse_error = 0;
    prop_error = 0;
    for numIter = 1:1:NUM_ITER
        
        beta_sum = 0;
        sum_G = zeros(M(m_idx),N);
        Gil = zeros(M(m_idx),K,L);
        for l=1:1:L
            
            % Generate channels.
            beta = a;
            if(l == 1)
                beta = 1;
            end
            betaMatrix = sqrt(beta)*eye(K);
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(M(m_idx),K),randn(M(m_idx),K))*betaMatrix;
            
            % Summation of all channels.
            sum_G = sum_G + Gil(:,:,l)*(S');
            
            % Summation of betas.
            beta_sum = beta_sum + beta;
            
        end
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q*N));
        
        % Apply squared pilot power.
        sum_G = sqrt(q)*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M(m_idx),N),randn(M(m_idx),N));
        
        % received pilot symbols at BS.
        Y1 = sum_G + W1;
        
        % ******* LS ********
        Z11_ls = (1/(sqrt(q)*N))*Y1*S(:,1);       
        ls_error = ls_error + ((Z11_ls'-Gil(:,1,1)')*(Z11_ls-Gil(:,1,1)));
        
        % ******* MMSE ********
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        mmse_error = mmse_error + ((Z11_mmse'-Gil(:,1,1)')*(Z11_mmse-Gil(:,1,1)));
        
        % ******* Prop. ********
        Z11_prop = M(m_idx)*beta111*Z11_ls/((Z11_ls')*(Z11_ls));
        prop_error = prop_error + ((Z11_prop'-Gil(:,1,1)')*(Z11_prop-Gil(:,1,1)));
    end
    
    fprintf('M: %d\n',M(m_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M(m_idx) * NUM_ITER);
    ls_error_vec(m_idx) = ls_error;
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error(m_idx) = epsilon11 - beta111;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M(m_idx) * NUM_ITER);
    mmse_error_vec(m_idx) = mmse_error;
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error(m_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
    
    % Proposed estimator Simulation Error (Monte Carlo)
    prop_error = prop_error./(M(m_idx) * NUM_ITER);
    prop_error_vec(m_idx) = prop_error;   
    
    % Proposed estimator Theoretical error.
    theta_ik = calculateTheta_ikv2(beta111, epsilon11, M(m_idx));
    theoretical_proposed_error(m_idx) = ((M(m_idx)/(M(m_idx)-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
    
    % Proposed estimator Approximated error.
    theoretical_proposed_approx_error(m_idx) = (beta111*((((2-M(m_idx))*beta111)/((M(m_idx)-1)*epsilon11)) + 1));
end

fontSize = 10;

fdee_figure = figure;
semilogy(M,real(theoretical_mmse_error),'-b');
hold on;
semilogy(M,real(mmse_error_vec),'b*','MarkerSize',7);
semilogy(M,real(theoretical_ls_error),'-g','MarkerSize',7);
semilogy(M,real(ls_error_vec),'g*','MarkerSize',7);
semilogy(M,real(theoretical_proposed_error),'-r','MarkerSize',7);
semilogy(M,real(theoretical_proposed_approx_error),'rx','MarkerSize',7);
semilogy(M,real(prop_error_vec),'ro','MarkerSize',7);
hold off
grid on;
xlabel('M')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (approx.)', 'Prop. (sim)', 'Location','northwest');
strText = sprintf('q = %1.0f dB, a = %1.2f',q,a);
x1 = M(1)+(30/100)*M(1);
y1 = 1.5;
text(x1,y1,strText,'FontSize', fontSize)

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure2_MSE_versus_number_of_antennas_K%d_L%d_N%d_v0_%s.fig',K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure2_MSE_versus_number_of_antennas_K%d_L%d_N%d_v0_%s.mat',K,L,N,timeStamp);
save(fileName);
