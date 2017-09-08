clear all;close all;clc

%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v10_20160804T083250.mat');
%load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v13_20170319T153046.mat');
load('channel_estimation_mse_vs_tx_snr_M30_K10_P20_L7_N223_v13_20170320T033804.mat');

%SNR = -12:2:10;

SNR_idx_start = find(SNR==-12);
SNR_idx_end = find(SNR==20);

SNR_aux = SNR(SNR_idx_start):2:SNR(SNR_idx_end);

fontSize = 12;
lineWidth = 1;
markerSize = 8;

fdee_figure = figure;

if(1)
    semilogy(SNR_aux,real(prop_1_sim_error_vec(SNR_idx_start:SNR_idx_end)),'-k*','MarkerSize',markerSize, 'LineWidth', lineWidth);
    hold on;
    semilogy(SNR_aux,real(prop_2_sim_error_vec(SNR_idx_start:SNR_idx_end)),'-kv','MarkerSize',markerSize, 'LineWidth', lineWidth);
    semilogy(SNR_aux,real(prop_3_sim_error_vec(SNR_idx_start:SNR_idx_end)),'-ko','MarkerSize',markerSize, 'LineWidth', lineWidth);
    semilogy(SNR_aux,real(prop_4_sim_error_vec(SNR_idx_start:SNR_idx_end)),'-ks','MarkerSize',markerSize, 'LineWidth', lineWidth);
    semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r','MarkerSize',markerSize, 'LineWidth', lineWidth);
end

if(0)
    semilogy(SNR_aux,real(prop_1_ana_error_vec(SNR_idx_start:SNR_idx_end)),'-k*','MarkerSize',7);
    hold on;
    semilogy(SNR_aux,real(prop_2_ana_error_vec(SNR_idx_start:SNR_idx_end)),'-kv','MarkerSize',7);
    semilogy(SNR_aux,real(prop_3_ana_error_vec(SNR_idx_start:SNR_idx_end)),'-ko','MarkerSize',7);
    semilogy(SNR_aux,real(prop_4_ana_error_vec(SNR_idx_start:SNR_idx_end)),'-ks','MarkerSize',7);
    semilogy(SNR_aux,theoretical_mmse_error(SNR_idx_start:SNR_idx_end),'--r*');
end

hold off
grid on;
axis([-12 20 0.23 0.301])
xlabel('SNR [dB]')
ylabel('MSE')
leg1 = legend('Prop. avg. = 1 (sim.)','Prop. avg. = 5 (sim.)','Prop. avg. = 10 (sim.)','Prop. avg. = 20 (sim.)','MMSE (analytic.)');


set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])