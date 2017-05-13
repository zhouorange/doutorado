clear all;clc

load('Massive_MU_MIMO_M_100_K_2_4_6_8_10_flat_fading.mat');

LineWidth = 2;
FontSize = 14;

% Set up a figure for visualizing BER results.
figure;
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[0.9e-5 1e-2]);
xlabel('Eb/No (dB)','FontSize',FontSize,'FontName','Times'); ylabel('BER','FontSize',FontSize,'FontName','Times'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name', 'OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink','FontSize',FontSize,'FontName','Times');

for k_idx=1:1:length(K)
    
    for m_idx=1:1:length(M)
        
        % Loop over selected EbNo points.
        for idx = 1:length(snr)

            % Plot results
            semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'ro','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
            if(k_idx==1)
                semilogy(EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'ks','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
            end
            set(gca,'FontSize',FontSize)
            set(gca,'FontName','Times')
            legend('MRC', 'MFB');
            drawnow;
        end
        
        semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'r-','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
        if(k_idx==1)
            semilogy(EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'k-','LineWidth',LineWidth,'MarkerFaceColor','white','MarkerSize',8);
        end
              
    end
end
hold off;
