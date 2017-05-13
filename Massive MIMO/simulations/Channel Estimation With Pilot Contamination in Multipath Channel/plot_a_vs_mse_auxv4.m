clear all;close all;clc

num_of_points = 16;

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v14_20170322T003110.mat');% 1st Figure
theoretical_mmse_error_m30 = theoretical_mmse_error(1:num_of_points);
theoretical_ls_error_m30 = theoretical_ls_error(1:num_of_points);
prop_4_ana_error_vec_m30 = prop_4_ana_error_vec(1:num_of_points);

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

load('channel_estimation_mse_vs_tx_snr_M90_K10_P20_L7_N223_v14_20170322T101609.mat'); % 2nd Figure
theoretical_mmse_error_m90 = theoretical_mmse_error(1:num_of_points);
theoretical_ls_error_m90 = theoretical_ls_error(1:num_of_points);
prop_4_ana_error_vec_m90 = prop_4_ana_error_vec(1:num_of_points);

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

a = a(1:num_of_points);

fdee_figure = figure;

subplot(1,2,1)
loglog(a,theoretical_mmse_error_m30,'k-','MarkerSize',5);
hold on;
loglog(a,theoretical_ls_error_m30,'kv-','MarkerSize',5);
loglog(a,real(prop_4_ana_error_vec_m30),'k*-','MarkerSize',5);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)','Location','northwest');
%axis([1e-3 1 0.0063 5.6])
axis([1e-3 1 -Inf 4])

subplot(1,2,2)
loglog(a,theoretical_mmse_error_m90,'k-','MarkerSize',5);
hold on;
loglog(a,theoretical_ls_error_m90,'kv-','MarkerSize',5);
loglog(a,real(prop_4_ana_error_vec_m90),'k*-','MarkerSize',5);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)','Location','northwest');
%axis([1e-3 1 0.0063 5.6])
axis([1e-3 1 -Inf 4])
