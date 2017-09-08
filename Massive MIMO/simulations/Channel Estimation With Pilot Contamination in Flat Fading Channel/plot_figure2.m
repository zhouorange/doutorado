clear all;close all;clc

load('Figure2_MSE_versus_number_of_antennas_K10_L7_N10_v0_20170618T184816.mat');

fontSize = 10;

fdee_figure = figure;
semilogy(M,real(theoretical_mmse_error),'r-');
hold on;
semilogy(M,real(mmse_error_vec),'r-.s','MarkerSize',7);
semilogy(M,real(theoretical_ls_error),'b-','MarkerSize',7);
semilogy(M,real(ls_error_vec),'b^-.','MarkerSize',7);
semilogy(M,real(theoretical_proposed_approx_error),'k-','MarkerSize',7);
semilogy(M,real(prop_error_vec),'ko-.','MarkerSize',7);
hold off
grid on;
xlabel('M')
ylabel('MSE')
legend('MMSE (analytical)','MMSE (simulated)','LS (analytical)', 'LS (simulated)', 'Prop. (approximated)', 'Prop. (simulated)', 'Location','northwest');
strText = sprintf('q = %1.0f dB, a = %1.2f',q,a);
x1 = 130;
y1 = 0.27;
text(x1,y1,strText,'FontSize', fontSize)
axis([M(1) M(end) 0.23 0.33])