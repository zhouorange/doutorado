clear all;close all;clc

rng(1)

q_db = -10; %-50:4:50;     % Uplink pilot power.
q = 10.^(q_db./10);

M = 10:10:100;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.

a = 0.05;           % Constant beta value.
beta111 = 1;

NUM_ITER = 10000;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

theoretical_proposed_error = zeros(1,length(q_db));
theoretical_proposed_approx_error = zeros(1,length(q_db));
error_mse = zeros(1,length(q_db));
for q_idx=1:1:length(q_db)
    
    beta = generateStaticBeta(a, beta111, L, K);
    
    % Summation of betas.
    [beta_sum, beta11k] = getBetaSum(beta, 1, L);
    
    % Factor.
    epsilon11 = (beta_sum + 1/(q(q_idx)*N));
    
    fprintf('q: %d\n',q_db(q_idx));
    
    % Proposed estimator Theoretical error.
    theta_ik = calculateTheta_ikv2(beta111, epsilon11, M);
    theoretical_proposed_error(q_idx) = ((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
    
    % Proposed estimator Approximated error.
    theoretical_proposed_approx_error(q_idx) = (beta111*((((2-M)*beta111)/((M-1)*epsilon11)) + 1));
    
    error_mse(q_idx) = abs(theoretical_proposed_error(q_idx) - theoretical_proposed_approx_error(q_idx));
end

fontSize = 10;

fdee_figure = figure;
semilogy(q_db,real(error_mse),'--k','MarkerSize',7);
grid on;
xlabel('SNR [dB]')
ylabel('MSE')
%legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (sim)', 'Prop. (approx.)');
%axis([SNR(1) SNR(length(SNR)) 0.18 5.1])
strText = sprintf('M = %d, a = %1.2f',M,a);
x1 = q_db(length(q_db))-12;
y1 = 1;
text(x1,y1,strText,'FontSize', fontSize)

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure5_comparison_closed_form_and_approximated_M%d_K%d_L%d_N%d_v0_%s.fig',M,K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure5_comparison_closed_form_and_approximated_M%d_K%d_L%d_N%d_v0_%s.mat',M,K,L,N,timeStamp);
save(fileName);