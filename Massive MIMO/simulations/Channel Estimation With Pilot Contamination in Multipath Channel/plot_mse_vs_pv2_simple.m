clear all;close all;clc

%load('channel_estimation_mse_vs_P_K10_M30_L7_N313_v2_20160805T201921.mat'); % SNR = -10
%load('channel_estimation_mse_vs_P_K10_M30_L7_N313_v2_20160805T194757.mat'); % SNR = 0

%load('channel_estimation_mse_vs_P_K10_M30_L7_N313_v2_20160805T212157.mat'); %SNR = 0

%load('channel_estimation_mse_vs_P_K10_M30_L7_N563_v2_20160806T085039.mat'); %SNR = 0

%load('channel_estimation_mse_vs_P_K10_M30_L7_N563_v2_20160807T104153.mat'); %SNR = 20

%load('channel_estimation_mse_vs_P_K10_M30_L7_N563_v3_20170319T104041.mat'); % SNR = 20

load('channel_estimation_mse_vs_P_K10_M30_L7_N563_v4_20170413T225310.mat'); %SNR = 0 dB

p_idx = length(P);
% Plot results.
yyaxis left
semilogy(P(1:p_idx),real(theoretical_mmse_error(1:p_idx)),'r-', ...
    P(1:p_idx),real(theoretical_ls_error(1:p_idx)),'ko-', 'MarkerSize',6);
xlabel('P'); 
ylabel('MSE');

grid on;

axis([P(1) P(length(P)) 0.231 0.377])

yyaxis right
plot(P(1:p_idx), N_vector(1:p_idx), 'b*-','MarkerSize',7);
ylabel('N');
legend('MMSE (ideal)', 'LS (ana.)','# Symbol Pilots');

set(gca,'XTickMode','manual');
set(gca,'XTick',[P(1):3:P(length(P))]);

fprintf(1,'SNR: %d dB\n', SNR);
fprintf(1,'M: %d antennas\n', M);