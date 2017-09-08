clear all;close all;clc

%load('Figure1_MSE_versus_TX_SNR_M100_K10_L7_N10_v0_20170622T211139.mat');
load('Figure1_MSE_versus_TX_SNR_M70_K10_L7_N10_v0_20170618T175056.mat')

fontSize = 10;

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'r-');
hold on;
semilogy(SNR,real(mmse_error_vec),'r-.s','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'b-','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b^-.','MarkerSize',7);
semilogy(SNR,real(theoretical_proposed_approx_error),'k-','MarkerSize',7);
semilogy(SNR,real(prop_error_vec),'k*-.','MarkerSize',7);
hold off
grid on;
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (analytical)','MMSE (simulated)','LS (analytical)', 'LS (simulated)', 'Prop. (approximated)', 'Prop. (simulated)');
axis([SNR(1) SNR(length(SNR)) 0.18 5.1])
strText = sprintf('M = %d, a = %1.2f',M,a);
x1 = SNR(length(SNR))-12;
y1 = 1.4;
text(x1,y1,strText,'FontSize', fontSize)
