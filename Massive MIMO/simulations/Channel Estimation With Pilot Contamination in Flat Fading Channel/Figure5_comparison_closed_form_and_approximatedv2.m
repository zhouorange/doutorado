clear all;close all;clc

rng(1)

q_db = -50:4:50;     % Uplink pilot power.
q = 10.^(q_db./10);

M = [10 50 100 200 500];             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.

a = 0.5;           % Constant beta value.
beta111 = 1;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

error_mse = zeros(length(M),length(q));
for m_idx=1:1:length(M)
    
    theoretical_proposed_error = zeros(1,length(q));
    theoretical_proposed_approx_error = zeros(1,length(q));
    for q_idx=1:1:length(q)
        
        beta = generateStaticBeta(a, beta111, L, K);
        
        % Summation of betas.
        [beta_sum, beta11k] = getBetaSum(beta, 1, L);
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        % Proposed estimator Theoretical error.
        theta_ik = calculateTheta_ikv2(beta111, epsilon11, M(m_idx));
        theoretical_proposed_error(q_idx) = ((M(m_idx)/(M(m_idx)-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
        
        % Proposed estimator Approximated error.
        theoretical_proposed_approx_error(q_idx) = (beta111*((((2-M(m_idx))*beta111)/((M(m_idx)-1)*epsilon11)) + 1));
        
        error_mse(m_idx,q_idx) = abs(theoretical_proposed_error(q_idx) - theoretical_proposed_approx_error(q_idx));
        
        fprintf('M: %d - q: %d\n',M(m_idx),q_db(q_idx));
    end
end

fontSize = 10;
markerSize = 7;

fdee_figure = figure;
semilogy(q_db,real(error_mse(1,:)),'MarkerSize',markerSize);
hold on
semilogy(q_db,real(error_mse(2,:)),'MarkerSize',markerSize);
semilogy(q_db,real(error_mse(3,:)),'MarkerSize',markerSize);
semilogy(q_db,real(error_mse(4,:)),'MarkerSize',markerSize);
semilogy(q_db,real(error_mse(5,:)),'MarkerSize',markerSize);
hold off
grid on;
xlabel('q [dB]')
ylabel('Distance between analytical and approximated errors')
legend('M = 10','M = 50','M = 100', 'M = 200', 'M = 500','Location','southeast');
axis([q_db(1) q_db(length(q_db)) 3e-10 0.003])

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure5_comparison_closed_form_and_approximated_K%d_L%d_N%d_v0_%s.fig',K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure5_comparison_closed_form_and_approximated_K%d_L%d_N%d_v0_%s.mat',K,L,N,timeStamp);
save(fileName);