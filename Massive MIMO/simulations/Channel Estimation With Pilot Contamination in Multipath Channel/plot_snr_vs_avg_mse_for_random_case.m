clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v25_20170403T133522.mat')
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v27_20170404T044038.mat');

startSNR = 1;
endSNR = length(SNR)-1;

fdee_figure = figure;
semilogy(SNR(startSNR:endSNR),theoretical_mmse_error_vec(startSNR:endSNR),'r');
hold on;
semilogy(SNR(startSNR:endSNR),theoretical_ls_error_vec(startSNR:endSNR),'b','MarkerSize',7);
semilogy(SNR(startSNR:endSNR),theoretical_prop_error_vec(startSNR:endSNR),'ko-','MarkerSize',7);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec(startSNR:endSNR)),'k*','MarkerSize',7);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec_hat1(startSNR:endSNR)),'k-.','MarkerSize',7,'LineWidth',0.8);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec_hat2(startSNR:endSNR)),'k--','MarkerSize',7);
hold off
grid on;
axis([-22 0 0.001 1]);
xlabel('SNR [dB]')
ylabel('avg. MSE')
legend('MMSE (ana)','LS (ana)','Prop. (ana)','Prop. (sim)','Prop. (sim) \sigma^{2}=0.01','Prop. (sim) \sigma^{2}=0.001');
