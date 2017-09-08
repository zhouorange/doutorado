clear all;close all;clc

load('flat_channel_estimation_mse_vs_tx_snr_M30_K10_L7_N10_v7_20170503T200941.mat');

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
axis([SNR(1) SNR(length(SNR)) 0.2 10.4])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (sim)', 'Prop. (approx.)');
