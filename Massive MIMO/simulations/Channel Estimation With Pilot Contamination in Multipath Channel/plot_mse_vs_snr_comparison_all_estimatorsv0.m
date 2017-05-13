clear all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170318T141307.mat');
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170319T115554.mat');

%SNR = -12:2:10;

%SNR_idx_start = find(SNR==-12);
%SNR_idx_end = find(SNR==20);

SNR_idx_start = find(SNR==-12);
SNR_idx_end = find(SNR==6); 

SNR_aux = SNR(SNR_idx_start):2:SNR(SNR_idx_end);

fdee_figure = figure;
semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r');
hold on;
semilogy(SNR_aux,real(mmse_error_vec(SNR_idx_start:SNR_idx_end)),'r*','MarkerSize',7);
semilogy(SNR_aux,theoretical_ls_error(SNR_idx_start:SNR_idx_end),'--b','MarkerSize',7);
semilogy(SNR_aux,real(ls_error_vec(SNR_idx_start:SNR_idx_end)),'b*','MarkerSize',7);
semilogy(SNR_aux,real(theoretical_prop_error(SNR_idx_start:SNR_idx_end)),'--k','MarkerSize',7);
semilogy(SNR_aux,real(prop_error_vec(SNR_idx_start:SNR_idx_end)),'-k*','MarkerSize',7);
hold off
grid on;
%axis([-12 20 0.23 0.3712])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)','Prop. (ana)', 'Prop. (sim)');