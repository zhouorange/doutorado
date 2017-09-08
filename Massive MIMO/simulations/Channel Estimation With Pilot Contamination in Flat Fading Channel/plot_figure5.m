clear all;close all;clc

load('Figure5_comparison_closed_form_and_approximated_K10_L7_N10_v0_20170618T202020.mat');

fontSize = 10;
markerSize = 7;

fdee_figure = figure;
semilogy(q_db,real(error_mse(1,:)),'^-','MarkerSize',markerSize);
hold on
semilogy(q_db,real(error_mse(2,:)),'o-','MarkerSize',markerSize);
semilogy(q_db,real(error_mse(3,:)),'x-','MarkerSize',markerSize);
semilogy(q_db,real(error_mse(4,:)),'s-','MarkerSize',markerSize);
semilogy(q_db,real(error_mse(5,:)),'+-','MarkerSize',markerSize);
hold off
grid on;
xlabel('SNR [dB]')
ylabel('Distance between analytical and approximated errors')
legend('M = 10','M = 50','M = 100', 'M = 200', 'M = 500','Location','southeast');
axis([q_db(1) q_db(length(q_db)) 3e-10 1e-2])
