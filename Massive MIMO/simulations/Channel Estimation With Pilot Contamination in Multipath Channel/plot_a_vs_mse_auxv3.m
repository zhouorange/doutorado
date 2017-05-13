clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170322T003110.mat');% 1st Figure
load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T101609.mat'); % 2nd Figure
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170322T102534.mat');
%load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T131549.mat');
%load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T184816.mat');

num_of_points = 16;

a = a(1:num_of_points);
theoretical_mmse_error = theoretical_mmse_error(1:num_of_points);
theoretical_ls_error = theoretical_ls_error(1:num_of_points);
prop_4_ana_error_vec = prop_4_ana_error_vec(1:num_of_points);

fdee_figure = figure;
loglog(a,theoretical_mmse_error,'--r');
hold on;
loglog(a,theoretical_ls_error,'--b','MarkerSize',7);
loglog(a,real(prop_4_ana_error_vec),'--c','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d - SNR: %d [dB]',M,K,P,L,N,SNR);
title(titleStr);
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)','Location','northwest');
%axis([1e-3 1 0.0063 5.6])
axis([1e-3 1 -Inf 4.5])
