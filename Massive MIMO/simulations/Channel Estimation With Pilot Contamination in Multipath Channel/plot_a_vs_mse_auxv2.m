clear all;clc

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170322T003110.mat');
%load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T101609.mat');
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170322T102534.mat');
%load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T131549.mat');
%load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T184816.mat');

fdee_figure = figure;
loglog(a,theoretical_mmse_error,'--r');
hold on;
loglog(a,real(mmse_error_vec),'r*','MarkerSize',7);
loglog(a,theoretical_ls_error,'--b','MarkerSize',7);
loglog(a,real(ls_error_vec),'b*','MarkerSize',7);
loglog(a,real(prop_1_ana_error_vec),'--k','MarkerSize',7);
loglog(a,real(prop_1_sim_error_vec),'k*','MarkerSize',7);
loglog(a,real(prop_2_ana_error_vec),'--g','MarkerSize',7);
loglog(a,real(prop_2_sim_error_vec),'g*','MarkerSize',7);
loglog(a,real(prop_3_ana_error_vec),'--m','MarkerSize',7);
loglog(a,real(prop_3_sim_error_vec),'m*','MarkerSize',7);
loglog(a,real(prop_4_ana_error_vec),'--c','MarkerSize',7);
loglog(a,real(prop_4_sim_error_vec),'c*','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d - SNR: %d [dB]',M,K,P,L,N,SNR);
title(titleStr);
xlabel('a')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)','Prop. 1 avg=1 (ana)','Prop. 1 avg=1 (sim)','Prop. 2 avg=5 (ana)','Prop. 2 avg=5 (sim)','Prop. 3 avg=10 (ana)','Prop. 3 avg=10 (sim)','Prop. 4 avg=20 (ana)','Prop. 4 avg=20 (sim)');
