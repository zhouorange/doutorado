clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170318T141307.mat');
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170319T115554.mat');

load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170507T170825.mat');

snrStart = -30;
snrEnd = 10;

SNR_idx_start = find(SNR==snrStart);
SNR_idx_end = find(SNR==snrEnd);

SNR_aux = SNR(SNR_idx_start):4:SNR(SNR_idx_end);

fontSize = 12;
lineWidth = 1;
markerSize = 8;

fdee_figure = figure;

subplot(2,1,1)
semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r','MarkerSize',markerSize, 'LineWidth', lineWidth);
hold on;
semilogy(SNR_aux,real(mmse_error_vec(SNR_idx_start:SNR_idx_end)),'r*','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,theoretical_ls_error(SNR_idx_start:SNR_idx_end),'--b','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,real(ls_error_vec(SNR_idx_start:SNR_idx_end)),'b*','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,real(theoretical_prop_error(SNR_idx_start:SNR_idx_end)),'--k','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,real(prop_error_vec(SNR_idx_start:SNR_idx_end)),'ko','MarkerSize',markerSize, 'LineWidth', lineWidth);
hold off
grid on;
axis([snrStart snrEnd 0.2 5])
xlabel('SNR [dB]')
ylabel('MSE')
leg1 = legend('MMSE (analytic.)','MMSE (sim.)','LS (analytic.)', 'LS (sim.)','Prop. (analytic.)', 'Prop. (sim.)');
title('(\it{a})')

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')


load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v12_20170507T203152.mat');

snrStart = -15;
snrEnd = 10;

SNR_idx_start = find(SNR==snrStart);
SNR_idx_end = find(SNR==snrEnd);

SNR_aux = SNR(SNR_idx_start):1:SNR(SNR_idx_end);

subplot(2,1,2)
semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r','MarkerSize',markerSize, 'LineWidth', lineWidth);
hold on;
semilogy(SNR_aux,real(mmse_error_vec(SNR_idx_start:SNR_idx_end)),'r*','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,real(theoretical_prop_error(SNR_idx_start:SNR_idx_end)),'--k','MarkerSize',markerSize, 'LineWidth', lineWidth);
semilogy(SNR_aux,real(prop_error_vec(SNR_idx_start:SNR_idx_end)),'k*','MarkerSize',markerSize, 'LineWidth', lineWidth);
hold off
grid on;
axis([-15 10 0.228 0.31])
xlabel('SNR [dB]')
ylabel('MSE')
leg2 = legend('MMSE (analytic.)','MMSE (sim.)','Prop. (analytic.)', 'Prop. (sim.)');
title('(\it{b})')

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg2, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])
