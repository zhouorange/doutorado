clear all;close all;clc

SNR = 10;                           % Signal-to-noise ratio in dB.
K = 10;                             % Number of single-antenna users.
P = 20;                             % Channel Length (Finite Impulse Response - FIR).
L = 7;                              % Number of cells.
a = logspace(-3, 1, 20);

num_of_points = 16;

M1 = 30;                            % Number of antennas.
[theoretical_mmse_error, theoretical_ls_error, prop_4_ana_error_vec] = function_mse_vs_a(SNR, M1, K, P, L, a);
theoretical_mmse_error_m30 = theoretical_mmse_error(1:num_of_points);
theoretical_ls_error_m30 = theoretical_ls_error(1:num_of_points);
prop_4_ana_error_vec_m30 = prop_4_ana_error_vec(1:num_of_points);

M2 = 100;                            % Number of antennas.
[theoretical_mmse_error, theoretical_ls_error, prop_4_ana_error_vec] = function_mse_vs_a(SNR, M2, K, P, L, a);
theoretical_mmse_error_m90 = theoretical_mmse_error(1:num_of_points);
theoretical_ls_error_m90 = theoretical_ls_error(1:num_of_points);
prop_4_ana_error_vec_m90 = prop_4_ana_error_vec(1:num_of_points);

a = a(1:num_of_points);

fontSize = 12;
lineWidth = 1;
markerSize = 8;

fdee_figure = figure;

subplot(1,2,1)
loglog(a,theoretical_mmse_error_m30,'k-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on;
loglog(a,theoretical_ls_error_m30,'bv-','MarkerSize', markerSize, 'LineWidth', lineWidth);
loglog(a,real(prop_4_ana_error_vec_m30),'r*-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold off
grid on;
xlabel('Cross-Cell Interference, \it{a}')
ylabel('MSE')
leg1 = legend('MMSE (analytic.)','LS (analytic.)','Prop. (analytic.)','Location','northwest');
txt1 = sprintf('M = %d, q = %d dB',M1,SNR);
x1 = 0.002;
y1 = 0.9;
tx1 = text(x1,y1,txt1,'FontSize', fontSize);
axis([1e-3 1 -Inf 4])

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')
set(tx1, 'FontName', 'Times New Roman')

subplot(1,2,2)
loglog(a,theoretical_mmse_error_m90,'k-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on;
loglog(a,theoretical_ls_error_m90,'bv-','MarkerSize', markerSize, 'LineWidth', lineWidth);
loglog(a,real(prop_4_ana_error_vec_m90),'r*-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold off
grid on;
xlabel('Cross-Cell Interference, \it{a}')
ylabel('MSE')
leg2 = legend('MMSE (analytic.)','LS (analytic.)','Prop. (analytic.)','Location','northwest');
txt1 = sprintf('M = %d, q = %d dB',M2,SNR);
x1 = 0.002;
y1 = 0.9;
tx2 = text(x1,y1,txt1,'FontSize', fontSize);
axis([1e-3 1 -Inf 4])

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg2, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')
set(tx2, 'FontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])
