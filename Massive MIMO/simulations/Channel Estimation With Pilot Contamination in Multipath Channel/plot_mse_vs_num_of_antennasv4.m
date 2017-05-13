clear all;close all;clc

load('channel_estimation_mse_vs_num_antennas_K10_P20_L7_N223_v2_20160802T204318.mat')

%m_idx = length(M)-7

% Set up a figure for visualizing BER results.
fdee_figure = figure; 
grid on; hold on;
xlabel('M'); 
ylabel('MSE');
set(gca,'yscale','log');

% Plot results.
semilogy(M(1:m_idx),real(mmse_error_vec(1:m_idx)),'r*', ...
    M(1:m_idx),real(ls_error_vec(1:m_idx)),'b*', ...
    M(1:m_idx),real(prop_1_error_vec(1:m_idx)),'ko', ...
    M(1:m_idx),real(prop_3_error_vec(1:m_idx)),'ks', ...
    M(1:m_idx),real(prop_4_error_vec(1:m_idx)),'k*' ,'MarkerSize',7)
legend('MMSE (ideal)', 'LS (ana)', 'Prop. avg=1 (sim)', 'Prop. avg=10 (sim)', 'Prop. avg=20 (sim)');

semilogy(M,real(mmse_error_vec),'r-', ...
    M,real(ls_error_vec),'b-', ...
    M,real(prop_1_error_vec),'k-', ...
    M,real(prop_3_error_vec),'k-', ...
    M,real(prop_4_error_vec),'k-' ,'MarkerSize',7)
hold off;

axis([M(1) M(end) 0.23 0.3615])