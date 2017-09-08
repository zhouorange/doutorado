clear all;close all;clc

load('Figure4_MSE_under_random_large_scale_fading_M30_K10_L7_N10_v0_20170517T231803.mat');

fontSize = 10;

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'b');
hold on;
semilogy(SNR,theoretical_ls_error,'^r-','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_approx_error),'ok-','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'xk','MarkerSize',7);
semilogy(SNR,real(prop_error_hat1_vec),'k:','MarkerSize',7);
semilogy(SNR,real(prop_error_hat2_vec),'k-.','MarkerSize',7);
hold off
grid on;
xlabel('SNR [dB]')
ylabel('avg. MSE')
legend('MMSE (analytical)','LS (analytical)','Prop. (approximated)', 'Prop. (simulated)', 'Prop. (simulated) \sigma^{2} = 0.1','Prop. (simulated) \sigma^{2} = 0.01');
axis([-10 7 0.001 1])
strText = sprintf('M = %d',M);
x1 = 2;
y1 = 0.07;
text(x1,y1,strText,'FontSize', fontSize)