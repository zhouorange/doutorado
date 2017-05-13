clear all;close all;clc

LineWidth = 2;
FontSize = 14;

marker_size = 8;

load('./MSE_Comparison_LS_MMSE_OMP.mat');

EbNoVec = EbNoVec_minus10_upto_100;

est_error_ls = est_error_ls_minus10_upto_100;
est_error_mmse = est_error_mmse_minus10_upto_100;
est_error_omp = est_error_omp_minus10_upto_100;
est_error_omp_Np17 = est_error_omp_minus10_upto_100_Np17;

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-9 100]);
xlabel('Eb/No (dB)'); ylabel('MSE'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); 
set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
%strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
%title(strTitle);

semilogy(EbNoVec, est_error_ls(1,:), 'ko-','LineWidth',2,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec, est_error_mmse(1,:), 'kx','LineWidth',2,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec, est_error_omp(1,:), 'ks','LineWidth',2,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec, est_error_omp_Np17(1,:),'k-','LineWidth',2,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');

set(gca,'FontSize',FontSize)
set(gca,'FontName','Times New Roman')
legend('LS,        N_p = 32 pilots','MMSE, N_p = 32 pilots','OMP,    N_p = 32 pilots', 'OMP,    N_p = 17 pilots');
hold off

% semilogy(EbNoVec, avg_error_ls(1,:), 'bs');
% semilogy(EbNoVec, avg_error_mmse(1,:), 'k*');
% semilogy(EbNoVec, avg_error_omp(1,:), 'ro');
% legend('LS','MMSE','OMP');
% hold off