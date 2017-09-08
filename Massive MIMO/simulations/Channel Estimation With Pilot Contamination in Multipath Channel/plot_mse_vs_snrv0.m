clear all;clc

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v10_20160804T083250.mat');

%SNR = -10:2:10;

SNR_idx_start = find(SNR==-10);
SNR_idx_end = find(SNR==20); 

SNR_aux = SNR(SNR_idx_start):2:SNR(SNR_idx_end);

fdee_figure = figure;
semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r');
hold on;
semilogy(SNR_aux,real(mmse_error_vec(SNR_idx_start:SNR_idx_end)),'r*','MarkerSize',7);
semilogy(SNR_aux,theoretical_ls_error(SNR_idx_start:SNR_idx_end),'--b','MarkerSize',7);
semilogy(SNR_aux,real(ls_error_vec(SNR_idx_start:SNR_idx_end)),'b*','MarkerSize',7);
semilogy(SNR_aux,theoretical_proposed_error(SNR_idx_start:SNR_idx_end),'--k','MarkerSize',7);
semilogy(SNR_aux,real(prop_1_error_vec(SNR_idx_start:SNR_idx_end)),'k*','MarkerSize',7);
semilogy(SNR_aux,real(prop_2_error_vec(SNR_idx_start:SNR_idx_end)),'-kv','MarkerSize',7);
semilogy(SNR_aux,real(prop_3_error_vec(SNR_idx_start:SNR_idx_end)),'-ko','MarkerSize',7);
semilogy(SNR_aux,real(prop_4_error_vec(SNR_idx_start:SNR_idx_end)),'-ks','MarkerSize',7);
hold off
grid on;
%axis([-10 22 0.2 1.304])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)','Prop. (ana)','Prop. 1 avg=1','Prop. 2 avg=5','Prop. 3 avg=10', 'Prop. 4 avg=20');