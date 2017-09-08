clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v10_20160804T083250.mat');
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v13_20170319T153046.mat');
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v13_20170320T033804.mat');

%SNR = -12:2:10;

SNR_idx_start = find(SNR==-12);
SNR_idx_end = find(SNR==20);

SNR_aux = SNR(SNR_idx_start):2:SNR(SNR_idx_end);

fdee_figure = figure;

if(1)
    semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'-r');
    hold on;
    semilogy(SNR_aux,mmse_error_vec(SNR_idx_start:SNR_idx_end),'r*');
    semilogy(SNR_aux,real(theoretical_ls_error(SNR_idx_start:SNR_idx_end)),'-k','MarkerSize',7);
    semilogy(SNR_aux,real(ls_error_vec(SNR_idx_start:SNR_idx_end)),'k*','MarkerSize',7);
end

hold off
grid on;
axis([-12 20 0.23 0.3715])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ideal)','MMSE (sim)','LS (ana)','LS (sim)');