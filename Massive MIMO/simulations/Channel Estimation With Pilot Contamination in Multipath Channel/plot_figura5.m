clear all;close all;clc

load('channel_estimation_mse_vs_num_antennas_K10_P20_L7_N223_v2_20160802T204318.mat')

fontSize = 12;
lineWidth = 1;
markerSize = 8;

% Set up a figure for visualizing BER results.
fdee_figure = figure; 
grid on; hold on;
xlabel('M'); 
ylabel('MSE');
set(gca,'yscale','log');

% Plot results.
semilogy(M(1:m_idx),real(mmse_error_vec(1:m_idx)),'r-', ...
    M(1:m_idx),real(ls_error_vec(1:m_idx)),'b^-', ...
    M(1:m_idx),real(prop_1_error_vec(1:m_idx)),'ko-', ...
    M(1:m_idx),real(prop_3_error_vec(1:m_idx)),'ks-', ...
    M(1:m_idx),real(prop_4_error_vec(1:m_idx)),'k*-' ,'MarkerSize', markerSize, 'LineWidth', lineWidth)
leg1 = legend('MMSE (analytic.)', 'LS (analytic.)', 'Prop. avg. = 1 (sim.)', 'Prop. avg. = 10 (sim.)', 'Prop. avg. = 20 (sim.)');

hold off;

axis([M(1) M(end) 0.23 0.3615])

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])