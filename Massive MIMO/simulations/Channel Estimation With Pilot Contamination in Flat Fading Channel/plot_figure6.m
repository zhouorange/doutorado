clear all;close all;clc

load('Figure6_remark_4_K10_L7_N10_v0_20170618T203748.mat');

semilogy(SNR,distance(1,:),'^r-');
hold on
semilogy(SNR,distance(2,:),'xb-');
semilogy(SNR,distance(3,:),'ok-');
semilogy(SNR,distance(4,:),'gs-');
semilogy(SNR,distance(5,:),'+c-');
xlabel('SNR[dB]')
ylabel('Distance between Prop. and MMSE estimators')
%axis([SNR(1) SNR(length(SNR)) 9e-6 0.005])
legend('M = 10','M = 30','M = 100','M = 200','M = 500','Location','southeast')
grid on
hold off
axis([SNR(1) 100 1e-4 1e-1])