clear all;close all;clc

%load('./Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150624T204448.mat');
%load('Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150625T080649.mat');
%load('Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150625T223301.mat');
%load('Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150625T233439.mat');
%load('Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150626T001840.mat');
%load('Massive_MU_MIMO_M_100_K_10_Np_101_estimation_MSE_20150626T003140.mat')

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('MMSE'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
title(strTitle);

semilogy(EbNoVec, est_error_ls(1,:), 'bs');
semilogy(EbNoVec, est_error_mmse(1,:), 'k*');
semilogy(EbNoVec, est_error_omp(1,:), 'ro');
legend('LS','MMSE','OMP');
hold off

% semilogy(EbNoVec, avg_error_ls(1,:), 'bs');
% semilogy(EbNoVec, avg_error_mmse(1,:), 'k*');
% semilogy(EbNoVec, avg_error_omp(1,:), 'ro');
% legend('LS','MMSE','OMP');
% hold off