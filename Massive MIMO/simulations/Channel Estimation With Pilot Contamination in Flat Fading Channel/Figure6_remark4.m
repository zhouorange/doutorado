clear all;close all;clc

SNR = -10:4:102;                           % Signal-to-Noise ratio, given in dB.
linear_SNR = 10.^(SNR./10);

M = 1e20;
%M = [10 30 100 200 500];         % Number of antennas.
K = 10;                             % Number of single-antenna users.
L = 7;                              % Number of cells.

N = K;            % Pilot length is set according to K and P.

a = 0.05;                           % Cross-cell interference level.
beta111 = 1;

beta = generateStaticBeta(a, beta111, L, K);

% Summation of betas.
[beta_sum, beta111] = getBetaSum(beta, 1, L);

sum_of_all_betas = getSumOfAllBetas(beta,L,K);

% ******** Simulation starts here. ********
distance = zeros(length(M),length(SNR));
for m_idx=1:1:length(M)
    for snr_idx = 1:1:length(SNR)
        
        q = (K*linear_SNR(snr_idx))./(N*(sum_of_all_betas));                % Uplink pilot power.
        
        epsilon_11 = (beta_sum + 1/(q*N));
        
        distance(m_idx,snr_idx) = (1./((M(m_idx))-1)).*(beta111.^2).*(1./epsilon_11);
    end
    fprintf(1, 'M: %d\n',M(m_idx));
end

fdee_figure = figure;
if(length(M) >= 1)
    semilogy(SNR,distance(1,:));
end
hold on
if(length(M) >= 2)
    semilogy(SNR,distance(2,:));
end
if(length(M) >= 3)
    semilogy(SNR,distance(3,:));
end
if(length(M) >= 4)
    semilogy(SNR,distance(4,:));
end
if(length(M) >= 5)
    semilogy(SNR,distance(5,:));
end
xlabel('SNR[dB]')
ylabel('Distance between Prop. and MMSE estimators')
%axis([SNR(1) SNR(length(SNR)) 9e-6 0.005])
legend('M = 10','M = 30','M = 100','M = 200','M = 500','Location','southeast')
grid on
hold off

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure6_remark_4_K%d_L%d_N%d_v0_%s.fig',K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure6_remark_4_K%d_L%d_N%d_v0_%s.mat',K,L,N,timeStamp);
save(fileName);
