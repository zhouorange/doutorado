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
    M(1:m_idx),real(ls_error_vec(1:m_idx)),'b*','MarkerSize',7)
legend('MMSE (ideal)', 'LS (ana)');

semilogy(M,real(mmse_error_vec),'r-', ...
    M,real(ls_error_vec),'b-','MarkerSize',7)
hold off;

axis([M(1) M(end) 0.23 0.3615])