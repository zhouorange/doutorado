clear all;close all;clc

load('Figure4_MSE_under_random_large_scale_fading_M30_K10_L7_N10_v0_20170517T231803.mat');

fontSize = 10;

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'--r');
hold on;
semilogy(SNR,real(mmse_error_vec),'r*','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'--b','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b*','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_error),'--k','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'k*','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_approx_error),'ko','MarkerSize',7);
semilogy(SNR,real(prop_error_hat1_vec),'k:','MarkerSize',7);
semilogy(SNR,real(prop_error_hat2_vec),'k-.','MarkerSize',7);
hold off
grid on;
xlabel('SNR [dB]')
ylabel('avg. MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (sim)', 'Prop. (approx.)','Prop. (sim) \sigma^{2} = 0.1','Prop. (sim) \sigma^{2} = 0.01');
axis([-10 7 0.003 0.5])
strText = sprintf('M = %d',M);
x1 = SNR(length(SNR))-12;
y1 = 1;
text(x1,y1,strText,'FontSize', fontSize)