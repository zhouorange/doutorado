clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v25_20170403T133522.mat')
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v27_20170404T044038.mat');

startSNR = 1;
endSNR = length(SNR)-1;

fontSize = 12;
lineWidth = 1;
markerSize = 8;

fdee_figure = figure;
semilogy(SNR(startSNR:endSNR),theoretical_mmse_error_vec(startSNR:endSNR),'r','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on;
semilogy(SNR(startSNR:endSNR),theoretical_ls_error_vec(startSNR:endSNR),'b','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR(startSNR:endSNR),theoretical_prop_error_vec(startSNR:endSNR),'ko-','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec(startSNR:endSNR)),'k*','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec_hat1(startSNR:endSNR)),'k-.','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR(startSNR:endSNR),real(prop_error_vec_hat2(startSNR:endSNR)),'k--','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold off
grid on;
axis([-22 0 0.001 1]);
xlabel('SNR [dB]')
ylabel('avg. MSE')
leg1 = legend('MMSE (analytic.)','LS (analytic.)','Prop. (analytic.)','Prop. (sim.)','Prop. (sim.) \sigma^{2} = 0.01','Prop. (sim.) \sigma^{2} = 0.001');

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])
