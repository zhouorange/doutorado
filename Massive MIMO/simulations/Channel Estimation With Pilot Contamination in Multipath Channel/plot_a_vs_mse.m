clear all;close all;clc

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170320T230618.mat');
theoretical_mmse_error_m30 = theoretical_mmse_error;
theoretical_ls_error_m30 = theoretical_ls_error;
prop_4_ana_error_vec_m30 = prop_4_ana_error_vec;

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170320T234116.mat');
theoretical_mmse_error_m90 = theoretical_mmse_error;
theoretical_ls_error_m90 = theoretical_ls_error;
prop_4_ana_error_vec_m90 = prop_4_ana_error_vec;

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

fdee_figure = figure;
subplot(1,2,1)
loglog(a,theoretical_ls_error_m30,'--b','MarkerSize',7);
hold on;
loglog(a,theoretical_mmse_error_m30,'--r');
loglog(a,real(prop_4_ana_error_vec_m30),'--c','MarkerSize',7);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('LS (ana)','MMSE (ideal)','Prop. (ana)');

subplot(1,2,2)
loglog(a,theoretical_ls_error_m90,'--b','MarkerSize',7);
hold on;
loglog(a,theoretical_mmse_error_m90,'--r');
loglog(a,real(prop_4_ana_error_vec_m90),'--c','MarkerSize',7);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('LS (ana)','MMSE (ideal)','Prop. (ana)');
