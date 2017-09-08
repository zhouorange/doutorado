clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v25_20170403T133522.mat')
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v27_20170404T044038.mat');
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v27_20170414T094819_simple.mat');

theoretical_mmse_error_vec1(1:18) = theoretical_mmse_error_vec;
mmse_error_vec1(1:18) = mmse_error_vec;
theoretical_ls_error_vec1(1:18) = theoretical_ls_error_vec;
ls_error_vec1(1:18) = ls_error_vec;
SNR1(1:18) = SNR;

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v27_20170415T004415_simple.mat');

total1 = length(theoretical_mmse_error_vec1);
total2 = length(theoretical_mmse_error_vec);
theoretical_mmse_error_vec1(total1+1:total1+total2) = theoretical_mmse_error_vec;
mmse_error_vec1(total1+1:total1+total2) = mmse_error_vec;
theoretical_ls_error_vec1(total1+1:total1+total2) = theoretical_ls_error_vec;
ls_error_vec1(total1+1:total1+total2) = ls_error_vec;
SNR1(total1+1:total1+total2) = SNR;

startSNR = 1;
endSNR = length(SNR1);

fdee_figure = figure;
semilogy(SNR1(startSNR:endSNR),theoretical_mmse_error_vec1(startSNR:endSNR),'r-');
hold on;
semilogy(SNR1(startSNR:endSNR),real(mmse_error_vec1(startSNR:endSNR)),'r*','MarkerSize',7);
semilogy(SNR1(startSNR:endSNR),theoretical_ls_error_vec1(startSNR:endSNR),'k-','MarkerSize',7);
semilogy(SNR1(startSNR:endSNR),real(ls_error_vec1(startSNR:endSNR)),'k*','MarkerSize',7);
hold off
grid on;
axis([-22 32 0.0002 0.712]);
xlabel('SNR [dB]')
ylabel('avg. MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)','LS (sim)');
