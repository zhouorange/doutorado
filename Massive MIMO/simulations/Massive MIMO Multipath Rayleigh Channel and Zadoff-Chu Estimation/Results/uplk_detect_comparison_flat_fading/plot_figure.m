


% Set up a figure for visualizing BER results.
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-6 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','MU-MIMO on Uplink');
title('OFDM modulated with QPSK System');

for idx = 1:length(EbNoVec)
    
    % Plot results
    semilogy(EbNoVec(1:idx), BER_ZF(1:idx), 'r*', ...
        EbNoVec(1:idx), BER_MMSE(1:idx), 'bo', ...
        EbNoVec(1:idx), BER_ML(1:idx), 'gs', ...
        EbNoVec(1:idx), BER_MRC(1:idx), 'ks', ...
        EbNoVec(1:idx), BER_ZF_LE(1:idx), 'b*', ...
        EbNoVec(1:idx), BER_MMSE_LE(1:idx), 'ko', ...
        EbNoVec(1:idx), BER_EGC(1:idx), 'go', ...
        EbNoVec(1:idx), BER_ZF_DF(1:idx), 'rs', ...
        EbNoVec(1:idx), BER_MMSE_DF(1:idx), 'bs');
    legend('ZF-SIC', 'MMSE-SIC', 'ML', 'MRC', 'ZF-LE', 'MMSE-LE', 'EGC', 'ZF-DF', 'MMSE-ZF');
    drawnow;
    
end

% Draw the lines
semilogy(EbNoVec, BER_ZF, 'r-', EbNoVec, BER_MMSE, 'b-', ...
    EbNoVec, BER_ML, 'g-', EbNoVec, BER_MRC, 'k-', EbNoVec, BER_ZF_LE, 'k-', EbNoVec, BER_MMSE_LE, 'k-', EbNoVec, BER_EGC, 'g-', EbNoVec, BER_ZF_DF, 'r-', EbNoVec, BER_MMSE_DF, 'b-');
hold off;






% Set up a figure for visualizing BER results.
figure;
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-6 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','MU-MIMO on Uplink');
title('OFDM modulated with QPSK System');

BER_MRC_smoothed = smooth(BER_MRC,'moving');
BER_EGC_smoothed = smooth(BER_EGC,'moving');

for idx = 1:length(EbNoVec)
    
    % Plot results
    semilogy(EbNoVec(1:idx), BER_ZF(1:idx), 'r*', ...
        EbNoVec(1:idx), BER_MMSE(1:idx), 'bo', ...
        EbNoVec(1:idx), BER_ML(1:idx), 'gs', ...
        EbNoVec(1:idx), BER_MRC_smoothed(1:idx), 'ks', ...
        EbNoVec(1:idx), BER_ZF_LE(1:idx), 'b*', ...
        EbNoVec(1:idx), BER_MMSE_LE(1:idx), 'ko', ...
        EbNoVec(1:idx), BER_EGC_smoothed(1:idx), 'go', ...
        EbNoVec(1:idx), BER_ZF_DF(1:idx), 'rs', ...
        EbNoVec(1:idx), BER_MMSE_DF(1:idx), 'bs');
    legend('ZF-SIC', 'MMSE-SIC', 'ML', 'MRC', 'ZF-LE', 'MMSE-LE', 'EGC', 'ZF-DF', 'MMSE-ZF');
    drawnow;
    
end

% Draw the lines
semilogy(EbNoVec, BER_ZF, 'r-', EbNoVec, BER_MMSE, 'b-', ...
    EbNoVec, BER_ML, 'g-', EbNoVec, BER_MRC_smoothed, 'k-', EbNoVec, BER_ZF_LE, 'k-', EbNoVec, BER_MMSE_LE, 'k-', EbNoVec, BER_EGC_smoothed, 'g-', EbNoVec, BER_ZF_DF, 'r-', EbNoVec, BER_MMSE_DF, 'b-');
hold off;