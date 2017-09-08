close all;clear all;clc

load('Figure3_MSE_versus_Cross_Cell_Interference_Level_M30_K10_L7_N10_v0_20170618T185955.mat');
theoretical_mmse_error_m30 = theoretical_mmse_error;
mmse_error_vec_3m0 = mmse_error_vec;
theoretical_ls_error_m30 = theoretical_ls_error;
ls_error_vec_m30 = ls_error_vec;
theoretical_proposed_approx_error_m30 = theoretical_proposed_approx_error;
prop_error_vec_m30 = prop_error_vec;
M1 = M;

load('Figure3_MSE_versus_Cross_Cell_Interference_Level_M90_K10_L7_N10_v0_20170618T184638.mat');
theoretical_mmse_error_m90 = theoretical_mmse_error;
mmse_error_vec_m90 = mmse_error_vec;
theoretical_ls_error_m90 = theoretical_ls_error;
ls_error_vec_m90 = ls_error_vec;
theoretical_proposed_approx_error_m90 = theoretical_proposed_approx_error;
prop_error_vec_m90 = prop_error_vec;
M2 = M;

fontSize = 10;
lineWidth = 1;
markerSize = 7;

fdee_figure = figure;

subplot(1,2,1)
loglog(a,theoretical_mmse_error_m30,'k-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on;
loglog(a,theoretical_ls_error_m30,'bv-','MarkerSize', markerSize, 'LineWidth', lineWidth);
loglog(a,real(theoretical_proposed_approx_error_m30),'r*-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold off
grid on;
xlabel('Cross-Cell Interference, \it{a}')
ylabel('MSE')
leg1 = legend('MMSE (analytical)','LS (analytical)','Prop. (approximated)','Location','northwest');
txt1 = sprintf('M = %d, q = %d dB',M1,SNR);
x1 = 0.015;
y1 = 1.6;
tx1 = text(x1,y1,txt1,'FontSize', fontSize);
axis([1e-2 1 -Inf 4])

% set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
% set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
% set(gca, 'defaultAxesFontName', 'Times New Roman')
% set(gca, 'defaultTextFontName', 'Times New Roman')
% set(tx1, 'FontName', 'Times New Roman')

subplot(1,2,2)
loglog(a,theoretical_mmse_error_m90,'k-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on;
loglog(a,theoretical_ls_error_m90,'bv-','MarkerSize', markerSize, 'LineWidth', lineWidth);
loglog(a,real(theoretical_proposed_approx_error_m90),'r*-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold off
grid on;
xlabel('Cross-Cell Interference, \it{a}')
ylabel('MSE')
leg2 = legend('MMSE (analytical)','LS (analytical)','Prop. (approximated)','Location','northwest');
txt1 = sprintf('M = %d, q = %d dB',M2,SNR);
x1 = 0.015;
y1 = 1.6;
tx2 = text(x1,y1,txt1,'FontSize', fontSize);
axis([1e-2 1 -Inf 4])

% set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
% set(leg2, 'FontName', 'Times New Roman', 'FontSize', fontSize)
% set(gca, 'defaultAxesFontName', 'Times New Roman')
% set(gca, 'defaultTextFontName', 'Times New Roman')
% set(tx2, 'FontName', 'Times New Roman')

% scaleFactor = 1.6;
% set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])
