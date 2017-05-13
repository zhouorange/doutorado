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

fdee_figure = figure;

subplot(1,2,1)
loglog(a,theoretical_mmse_error_m30,'k-','MarkerSize',5);
hold on;
loglog(a,theoretical_ls_error_m30,'bv-','MarkerSize',5);
loglog(a,real(prop_4_ana_error_vec_m30),'r*-','MarkerSize',5);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)','Location','northwest');
%axis([1e-3 1 0.0063 5.6])
txt1 = sprintf('M = %d, q = %d dB',M1,SNR);
%txt1 = 'M = 100, q = 10 dB';
x1 = 0.002;
y1 = 0.9;
text(x1,y1,txt1,'FontSize', 9)
axis([1e-3 1 -Inf 4])

subplot(1,2,2)
loglog(a,theoretical_mmse_error_m90,'k-','MarkerSize',5);
hold on;
loglog(a,theoretical_ls_error_m90,'bv-','MarkerSize',5);
loglog(a,real(prop_4_ana_error_vec_m90),'r*-','MarkerSize',5);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE (ideal)','LS (ana)','Prop. (ana)','Location','northwest');
txt1 = sprintf('M = %d, q = %d dB',M2,SNR);
%txt1 = 'M = 100, q = 10 dB';
x1 = 0.002;
y1 = 0.9;
text(x1,y1,txt1,'FontSize', 9)
%axis([1e-3 1 0.0063 5.6])
axis([1e-3 1 -Inf 4])
