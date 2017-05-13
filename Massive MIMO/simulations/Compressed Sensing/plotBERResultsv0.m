clear all;close all;clc

LineWidth = 1.5;
FontSize = 14;

marker_size = 8;

stepp = 1;

load('./BER_Comparison_LS_MMSE_OMP_M_100_K_10_Np_37.mat');

EbNoVec = EbNoVec_minus20_upto_40;

ber_mmse_ls = BER_MMSE_LS_minus20_upto_40;
ber_mmse_mmse = BER_MMSE_MMSE_minus20_upto_40;
ber_mmse_omp = BER_MMSE_OMP_minus20_upto_40;
ber_mmse_opt = BER_MMSE_OPT_minus20_upto_40;


% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[0.9e-5 0.5]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); 
set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
%strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
%title(strTitle);

semilogy(EbNoVec(1:stepp:end), ber_mmse_ls(1,1:stepp:end), 'ko-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec(1:stepp:end), ber_mmse_mmse(1,1:stepp:end), 'kx-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec(1:stepp:end), ber_mmse_omp(1,1:stepp:end), 'ks-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
semilogy(EbNoVec(1:stepp:end), ber_mmse_opt(1,1:stepp:end),'k-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');


% % Set up a figure for visualizing BER results.
% figura = figure; grid on; hold on;
% set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(9)],'ylim',[1e-5 0.5]);
% xlabel('Eb/No (dB)'); ylabel('BER'); set(figura,'NumberTitle','off');
% set(figura, 'renderer', 'zbuffer'); 
% set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
% %strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
% %title(strTitle);
% semilogy(EbNoVec(1:stepp:9), ber_mmse_ls(1,1:stepp:9), 'ko-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
% semilogy(EbNoVec(1:stepp:9), ber_mmse_mmse(1,1:stepp:9), 'kx-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
% semilogy(EbNoVec(1:stepp:9), ber_mmse_omp(1,1:stepp:9), 'ks-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');
% semilogy(EbNoVec(1:stepp:9), ber_mmse_opt(1,1:stepp:9),'k-','LineWidth',LineWidth,'MarkerSize',marker_size,'MarkerEdgeColor','k','MarkerFaceColor','w');

set(gca,'FontSize',FontSize)
set(gca,'FontName','Times New Roman')
legend('LS','MMSE','OMP', 'Full Channel');
hold off

% semilogy(EbNoVec, avg_error_ls(1,:), 'bs');
% semilogy(EbNoVec, avg_error_mmse(1,:), 'k*');
% semilogy(EbNoVec, avg_error_omp(1,:), 'ro');
% legend('LS','MMSE','OMP');
% hold off