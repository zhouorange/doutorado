clear all;clc

%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_100_K_10_flat_fading_various_detectors_20150208T141007.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_50_K_10_flat_fading_various_detectors_20150208T163418.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_250_K_10_flat_fading_various_detectors_20150208T222131.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_100_K_10_flat_fading_various_detectors_20150208T141007.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_50_K_10_flat_fading_various_detectors_20150208T163418.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_400_K_10_flat_fading_various_detectors_20150211T122801.mat');
%load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_500_K_10_flat_fading_various_detectors_20150212T174718.mat');
load('C:\Users\felipep\Dropbox\Doutorado\Massive MIMO\simulations\uplink_detection_comparison_flat_fading_results\Massive_MU_MIMO_M_400_K_10_flat_fading_various_detectors_20150211T122801.mat');

LineWidth = 2;
FontSize = 14;
LOC = 'SouthWest';

% Set up a figure for visualizing BER results.
figure;
h = gcf; grid on; hold on;
set(h, 'Position', [0 0 600 500])
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[0.9e-6 1e-1]);
xlabel('Eb/No (dB)','FontSize',FontSize,'FontName','Times'); ylabel('BER','FontSize',FontSize,'FontName','Times'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name', 'OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO Detection on Uplink','FontSize',FontSize,'FontName','Times');

for k_idx=1:1:length(K)
    
    numSym = K(k_idx)*NFFT; % number of symbols, i.e., number of terminals.
    
    for m_idx=1:1:length(M)
        
        % Loop over selected EbNo points.
        for idx = 1:length(EbNoVec)
            
            % Plot results
            semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'ro', EbNoVec(1:idx), BER_ZF_LE(m_idx,1:idx), 'bo', EbNoVec(1:idx), BER_MMSE_LE(m_idx,1:idx), 'ko', EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'ks','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
            set(gca,'FontSize',FontSize)
            set(gca,'FontName','Times')
            legend('MRC', 'ZF-LE', 'MMSE-LE', 'MFB','Location',LOC);
            drawnow;
        end
        
        semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'r-', EbNoVec(1:idx), BER_ZF_LE(m_idx,1:idx), 'b-', EbNoVec(1:idx), BER_MMSE_LE(m_idx,1:idx), 'k-', EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'k-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',8);
        
    end
end
hold off;

% Set up a figure for visualizing BER results.
figure;
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[-12.005, -11.93],'ylim',[2.1e-4 4.5e-4]);
xlabel('Eb/No (dB)','FontSize',FontSize,'FontName','Times'); ylabel('BER','FontSize',FontSize,'FontName','Times'); set(h,'NumberTitle','off');

for k_idx=1:1:length(K)
    
    numSym = K(k_idx)*NFFT; % number of symbols, i.e., number of terminals.
    
    for m_idx=1:1:length(M)
        
        % Loop over selected EbNo points.
        for idx = 1:length(snr)
            
            % Plot results
            semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'ro', EbNoVec(1:idx), BER_ZF_LE(m_idx,1:idx), 'bo', EbNoVec(1:idx), BER_MMSE_LE(m_idx,1:idx), 'ko', EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'ks','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
            set(gca,'FontSize',FontSize)
            set(gca,'FontName','Times')
            drawnow;
        end
        
        semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'r-', EbNoVec(1:idx), BER_ZF_LE(m_idx,1:idx), 'b-', EbNoVec(1:idx), BER_MMSE_LE(m_idx,1:idx), 'k-', EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'k-','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
        
    end
end
hold off;