clear all;clc

%---------- Parameters ----------
L = 1;                                                              % Number of cells.
K = 10;                                                             % Number of terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
modOrd = 2;                                                         % constellation size = 2^modOrd.
NFFT = 2048;                                                        % Number of points used by the OFDM.
NCP = 128;                                                          % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
EbNoVec = -25:1:0;                                                  % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

nTotalOfBits = 1e7;
nErrors = 100000;

%---------- Set up the simulation ----------
% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

% Preallocate variables for speed.
[BER_MRC, BER_MFB, BER_ZF_SIC, BER_MMSE_SIC, BER_ZF_LE, BER_MMSE_LE, BER_EGC, BER_ZF_DF, BER_MMSE_DF] = deal(zeros(length(M), length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[0.9e-7 1e-1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name', 'OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink');

for k_idx=1:1:length(K)
    
    numSym = K(k_idx)*NFFT; % number of symbols, i.e., number of terminals.
    
    for m_idx=1:1:length(M)
        
        % Loop over selected EbNo points.
        for idx = 1:length(snr)
            
            linearSnr = 10^(0.1*snr(idx));
            
            nErrs_mrc = 0;
            nErrs_mfb = 0;
            nErrs_zf_sic = 0;
            nErrs_mmse_sic = 0;
            nErrs_zf_le = 0;
            nErrs_mmse_le = 0;
            nErrs_egc = 0;
            nErrs_zf_df = 0;
            nErrs_mmse_df = 0;
            nBits = 0;
            nBits_mfb = 0;
            iter = 0;
            
            while(((nErrs_mrc < nErrors) || (nErrs_mfb < nErrors) || (nErrs_zf_sic < nErrors) || (nErrs_mmse_sic < nErrors) || (nErrs_zf_le < nErrors) || (nErrs_mmse_le < nErrors) || (nErrs_egc < nErrors) || (nErrs_zf_df < nErrors) || (nErrs_mmse_df < nErrors)) && (nBits_mfb < nTotalOfBits))
                
                iter = iter + 1;
                
                %% ---------- Transmission (UE) ----------
                % Create array of bits to modulate.
                msg = randi(hStr, [0 1], modOrd, numSym);
                
                msg_mfb = msg(:,1:K(k_idx):end);
                
                % Modulate data.
                source = modulate(hMod, msg);
                
                % Split source among K terminals.
                Tx = reshape(source, K(k_idx), numel(source)/K(k_idx)); clear source;
                
                % Create OFDM symbol.
                sequence = sqrt(NFFT)*ifft(Tx,NFFT,2);
                
                % Add CP.
                ofdm = [sequence(:,NFFT-NCP+1:end), sequence];
                
                % Make sure the OFDM symbol has unit variance.
                std_dev_vector = (sqrt(diag((ofdm*ofdm')/(NFFT+NCP))));
                std_dev_matrix = diag(1./std_dev_vector);
                ofdm = std_dev_matrix*ofdm;
                
                %% ---------- Channel ----------
                % Flat Rayleigh Fading - independent links.
                RayleighMat = (randn(hStr, M(m_idx), K(k_idx)) +  1i*randn(hStr, M(m_idx), K(k_idx)))/sqrt(2);
                
                % Add channel noise power to faded data.
                r = awgn(RayleighMat*ofdm, snr(idx), 0, hStr);
                
                % Add channel noise power to faded single user data.
                r_mfb = awgn(RayleighMat(:,1)*ofdm(1,:), snr(idx), 0, hStr);
                
                %% ---------- Reception (base station) ----------
                
                % Remove CP.
                rx = r(:,NCP+1:end);
                
                % Retrieve modulation symbols.
                ofdm_rx = (1/sqrt(NFFT))*fft(rx,NFFT,2);
                
                %% *********** 1) Maximum Ratio Combining (MRC) or Matched Filter *************
                E_mrc = zeros(modOrd, numSym);
                
                G = RayleighMat';
                r_mf = G*ofdm_rx;
                
                % Iterate over all Tx antennas and demodulate signal.
                for jj=1:1:K(k_idx)
                    demoded = demodulate(hDemod, r_mf(jj,:));
                    E_mrc(:,jj:K(k_idx):end) = reshape(demoded, modOrd, numSym/K(k_idx));
                end
                %******************************************************************************
                
                %% ********* 2) Matched Filter Bound (MFB) or Single User detection ***********
                
                % Remove CP.
                rx_mfb = r_mfb(:,NCP+1:end);
                
                % Retrieve modulation symbols.
                ofdm_rx_mfb = (1/sqrt(NFFT))*fft(rx_mfb,NFFT,2);
                
                % Apply MRC to direct path received signal.
                G = RayleighMat(:,1)';
                r_eq_mfb = G*ofdm_rx_mfb;
                
                % Demodulate signal.
                E_mfb = demodulate(hDemod, r_eq_mfb);
                % *****************************************************************************
                
                %% *********** 3) Zero-Forcing with Optimally Ordered SIC receiver ************
                E_zf_sic = zeros(modOrd, numSym); k = zeros(K(k_idx), 1);
                
                % Move OFDM symbol to aux variable.
                ofdm_rx_zf_sic = ofdm_rx;
                
                % Move the channel to aux variable.
                H = RayleighMat;
                
                % Initialization
                G = pinv(H);
                [val, k0] = min(sum(abs(G).^2,2));

                % Start Zero-Forcing Nulling Loop.
                for n = 1:K(k_idx)
                    % Find best transmitter signal using minimum norm.
                    k(n) = k0;
                    
                    % Select Weight vector for best transmitter signal.
                    w = G(k(n),:);
                    
                    % Calculate output for transmitter n.
                    y = w * ofdm_rx_zf_sic;
                    
                    % Demodulate bitstream.
                    demoded_zf_sic = demodulate(hDemod, y);
                    E_zf_sic(:, k(n):K(k_idx):end) = reshape(demoded_zf_sic, modOrd, numSym/K(k_idx));
                    
                    % Subtract effect of the transmitter n from received signal.
                    z = modulate(hMod, demodulate(hDemod, y));
                    ofdm_rx_zf_sic = ofdm_rx_zf_sic - H(:, k(n))*z;
                    
                    % Adjust channel estimate matrix for next minimum norm search.
                    H(:, k(n)) = zeros(M(m_idx), 1);
                    G = pinv(H);
                    for aa = 1:n
                        G(k(aa), :) = inf;
                    end
                    [val, k0] = min(sum(abs(G).^2,2));
                end
                % *****************************************************************************
                
                %% *************** 4) MMSE with Optimally Ordered SIC receiver ****************
                E_mmse_sic = zeros(modOrd, numSym); k = zeros(K(k_idx), 1);
                
                % Move OFDM symbol to aux variable.
                ofdm_rx_mmse_sic = ofdm_rx;
                
                % Move the channel to aux variable.
                H = RayleighMat;
                
                % Initialization
                G = (H*H' + (1/linearSnr)*eye(M(m_idx)))^-1;
                G = H'*G;
                [val, k0] = min(sum(abs(G).^2,2));
                % Start MMSE Nulling Loop
                for n = 1:K(k_idx)
                    % Find best transmitter signal using Min Norm
                    k(n) = k0;
                    
                    % Select Weight vector for best transmitter signal
                    w = G(k(n),:);
                    
                    % Calculate output for transmitter n and demodulate bitstream.
                    y = w * ofdm_rx_mmse_sic;
                    E_mmse_sic(:, k(n):K(k_idx):end) = reshape(demodulate(hDemod, y), modOrd, numSym/K(k_idx));
                    
                    % Subtract effect of the transmitter n from received signal
                    z = modulate(hMod, demodulate(hDemod, y));
                    ofdm_rx_mmse_sic = ofdm_rx_mmse_sic - H(:, k(n))*z;
                    
                    % Adjust channel estimate matrix for next min Norm search.
                    H(:, k(n)) = zeros(M(m_idx), 1);
                    G = (H*H' + (1/linearSnr)*eye(M(m_idx)))^-1;
                    G = H'*G;
                    for aa = 1:n
                        G(k(aa), :) = inf;
                    end
                    [val, k0] = min(sum(abs(G).^2,2));
                end
                % *****************************************************************************
                
                %% ***************************** 5) ZF-LE receiver ****************************
                % Initialization
                E_zf_le = zeros(modOrd, numSym);
                
                % Apply ZF to received symbol.
                G =  ((RayleighMat'*RayleighMat)^(-1))*RayleighMat';
                y = G*ofdm_rx;
                
                % Iterate over all Tx antennas and demodulate.
                for jj=1:1:K(k_idx)
                    demoded_zf_le = demodulate(hDemod, y(jj,:));
                    E_zf_le(:,jj:K(k_idx):end) = reshape(demoded_zf_le, modOrd, numSym/K(k_idx));
                end
                % *****************************************************************************
                
                %% ***************************** 6) MMSE-LE receiver ***************************
                % Initialization
                E_mmse_le = zeros(modOrd, numSym);
                
                % Move the channel to aux variable.
                H = RayleighMat;
                
                % Apply MMSE to received symbol.               
                G = (H*H' + (1/linearSnr)*eye(M(m_idx)))^-1;
                G = H'*G;
                y = G*ofdm_rx;
                
                % Iterate over all Tx antennas and demodulate.
                for jj=1:1:K(k_idx)
                    demoded_mmse_le = demodulate(hDemod, y(jj,:));
                    E_mmse_le(:,jj:K(k_idx):end) = reshape(demoded_mmse_le, modOrd, numSym/K(k_idx));
                end
                % *****************************************************************************
                
                %% ********************** 7) Equal Gain Combining (EGC) ************************
                % Initialization
                E_egc = zeros(modOrd, numSym);
                
                % Apply EGC to received OFDm signal.
                G = exp(-1*1i*angle(RayleighMat)).';
                
                % Removing the phase of the channel.
                yHat = G*ofdm_rx;
                
                % Iterate over all Tx antennas.
                for jj=1:1:K(k_idx)
                    demoded_egc = demodulate(hDemod, yHat(jj,:));
                    E_egc(:,jj:K(k_idx):end) = reshape(demoded_egc, modOrd, numSym/K(k_idx));
                end
                % *****************************************************************************
                
                %% ************* 8) Zero-Forcing Decision Feedback (ZF-DF) receiver ***********
                E_zf_df = zeros(modOrd, numSym);
                
                % Move OFDM symbol to aux variable.
                ofdm_rx_zf_df = ofdm_rx;
                
                % Move the channel to aux variable.
                H = RayleighMat;
                
                % Initialization
                G = pinv(H);
                % Start Zero-Forcing Nulling Loop.
                for n = 1:K(k_idx)
                    % Select Weight vector for best transmitter signal.
                    w = G(n,:);
                    
                    % Calculate output for transmitter n.
                    y = w * ofdm_rx_zf_df;
                    
                    % Demodulate bitstream.
                    demoded_zf_df = demodulate(hDemod, y);
                    E_zf_df(:, n:K(k_idx):end) = reshape(demoded_zf_df, modOrd, numSym/K(k_idx));
                    
                    % Subtract effect of the transmitter n from received signal.
                    z = modulate(hMod, demodulate(hDemod, y));
                    ofdm_rx_zf_df = ofdm_rx_zf_df - H(:, n)*z;
                end
                % *****************************************************************************
                               
                %% ************** 9) MMSE-Decision Feedback receiver (MMSE-DF) ****************
                E_mmse_df = zeros(modOrd, numSym);
                
                % Move OFDM symbol to aux variable.
                ofdm_rx_mmse_df = ofdm_rx;
                
                % Move the channel to aux variable.
                H = RayleighMat;                
                
                % Initialization.               
                G = (H*H' + (1/linearSnr)*eye(M(m_idx)))^-1;
                G = H'*G; 
                %G = (H'*H + ((1/linearSnr)*eye(K(k_idx)))) \ H';
                % Start MMSE Nulling Loop.
                for n = 1:K(k_idx)
                    % Select Weight vector for best transmitter signal.
                    w = G(n,:);
                    
                    % Calculate output for transmitter n and demodulate bitstream.
                    y = w * ofdm_rx_mmse_df;
                    E_mmse_df(:, n:K(k_idx):end) = reshape(demodulate(hDemod, y), modOrd, numSym/K(k_idx));
                    
                    % Subtract effect of the transmitter n from received signal.
                    z = modulate(hMod, demodulate(hDemod, y));
                    ofdm_rx_mmse_df = ofdm_rx_mmse_df - H(:, n)*z;
                end
                % *****************************************************************************
                
                %% ----------- Collecting Errors ---------------
                nErrs_mrc = nErrs_mrc + biterr(msg, E_mrc);
                nErrs_mfb = nErrs_mfb + biterr(msg_mfb, E_mfb);
                nErrs_zf_sic = nErrs_zf_sic + biterr(msg, E_zf_sic);
                nErrs_mmse_sic = nErrs_mmse_sic + biterr(msg, E_mmse_sic);
                nErrs_zf_le = nErrs_zf_le + biterr(msg, E_zf_le);
                nErrs_mmse_le = nErrs_mmse_le + biterr(msg, E_mmse_le);
                nErrs_egc = nErrs_egc + biterr(msg, E_egc);               
                nErrs_zf_df = nErrs_zf_df + biterr(msg, E_zf_df);
                nErrs_mmse_df = nErrs_mmse_df + biterr(msg, E_mmse_df);
                
                nBits = nBits + length(msg(:));
                nBits_mfb = nBits_mfb + length(msg_mfb(:));
                fprintf(1,'SNR: %d\n',EbNoVec(idx));
                fprintf(1,'BER_MRC: %f - nErrs_mrc: %d - nBits: %d - iter: %d\n',(nErrs_mrc./nBits),nErrs_mrc,nBits,iter);
                fprintf(1,'BER_MFB: %f - nErrs_mfb: %d - nBits_mfb: %d - iter: %d\n',(nErrs_mfb./nBits_mfb),nErrs_mfb,nBits_mfb,iter);
                fprintf(1,'BER_ZF_SIC: %f - nErrs_zf_sic: %d - nBits: %d - iter: %d\n',(nErrs_zf_sic./nBits),nErrs_zf_sic,nBits,iter);
                fprintf(1,'BER_MMSE_SIC: %f - nErrs_mmse_sic: %d - nBits: %d - iter: %d\n',(nErrs_mmse_sic./nBits),nErrs_mmse_sic,nBits,iter);
                fprintf(1,'BER_ZF_LE: %f - nErrs_zf_le: %d - nBits: %d - iter: %d\n',(nErrs_zf_le./nBits),nErrs_zf_le,nBits,iter);
                fprintf(1,'BER_MMSE_LE: %f - nErrs_mmse_le: %d - nBits: %d - iter: %d\n',(nErrs_mmse_le./nBits),nErrs_mmse_le,nBits,iter);
                fprintf(1,'BER_EGC: %f - nErrs_egc: %d - nBits: %d - iter: %d\n',(nErrs_egc./nBits),nErrs_egc,nBits,iter);               
                fprintf(1,'BER_ZF_DF: %f - nErrs_zf_df: %d - nBits: %d - iter: %d\n',(nErrs_zf_df./nBits),nErrs_zf_df,nBits,iter);
                fprintf(1,'BER_MMSE_DF: %f - nErrs_mmse_df: %d - nBits: %d - iter: %d\n',(nErrs_mmse_df./nBits),nErrs_mmse_df,nBits,iter);
                
                fprintf(1,'\n');
                
            end
            
            % Calculate BER for current point
            BER_MRC(m_idx,idx) = nErrs_mrc./nBits;
            BER_MFB(m_idx,idx) = nErrs_mfb./nBits_mfb;
            BER_ZF_SIC(m_idx,idx) = nErrs_zf_sic./nBits;
            BER_MMSE_SIC(m_idx,idx) = nErrs_mmse_sic./nBits;
            BER_ZF_LE(m_idx,idx) = nErrs_zf_le./nBits;
            BER_MMSE_LE(m_idx,idx) = nErrs_mmse_le./nBits;
            BER_EGC(m_idx,idx) = nErrs_egc./nBits;
            BER_ZF_DF(m_idx,idx) = nErrs_zf_df./nBits;
            BER_MMSE_DF(m_idx,idx) = nErrs_mmse_df./nBits;            
            
            % Plot results
            semilogy(EbNoVec(1:idx), BER_MRC(m_idx,1:idx), 'ro', EbNoVec(1:idx), BER_MFB(m_idx,1:idx), 'ks', EbNoVec(1:idx), BER_ZF_SIC(m_idx,1:idx), 'bs', EbNoVec(1:idx), BER_MMSE_SIC(m_idx,1:idx), 'rs', EbNoVec(1:idx), BER_ZF_LE(m_idx,1:idx), 'bo', EbNoVec(1:idx), BER_MMSE_LE(m_idx,1:idx), 'ko', EbNoVec(1:idx), BER_EGC(m_idx,1:idx), 'b*', EbNoVec(1:idx), BER_ZF_DF(m_idx,1:idx), 'r*', EbNoVec(1:idx), BER_MMSE_DF(m_idx,1:idx), 'k*');
            legend('MRC', 'MFB', 'ZF-SIC', 'MMSE-SIC', 'ZF-LE', 'MMSE-LE', 'EGC', 'ZF-DF', 'MMSE-DF');
            drawnow;
            
        end
        
        % Perform curve fitting and replot the results
        fitBER_MRC = berfit(EbNoVec, BER_MRC(m_idx,:));
        fitBER_MFB = berfit(EbNoVec, BER_MFB(m_idx,:));
        fitBER_ZF_SIC = berfit(EbNoVec, BER_ZF_SIC(m_idx,:));
        fitBER_MMSE_SIC = berfit(EbNoVec, BER_MMSE_SIC(m_idx,:));
        fitBER_ZF_LE = berfit(EbNoVec, BER_ZF_LE(m_idx,:));
        fitBER_MMSE_LE = berfit(EbNoVec, BER_MMSE_LE(m_idx,:));
        fitBER_EGC = berfit(EbNoVec, BER_EGC(m_idx,:));        
        fitBER_ZF_DF = berfit(EbNoVec, BER_ZF_DF(m_idx,:));
        fitBER_MMSE_DF = berfit(EbNoVec, BER_MMSE_DF(m_idx,:));       
        if(isempty(fitBER_MRC) || isempty(fitBER_MFB) || isempty(fitBER_ZF_SIC) || isempty(fitBER_MMSE_SIC) || isempty(fitBER_ZF_LE) || isempty(fitBER_MMSE_LE) || isempty(fitBER_EGC) || isempty(fitBER_ZF_DF) || isempty(fitBER_MMSE_DF) )
            semilogy(EbNoVec, BER_MRC(m_idx,:), 'r-', EbNoVec, BER_MFB(m_idx,:), 'k-', EbNoVec, BER_ZF_SIC(m_idx,:), 'b-', EbNoVec, BER_MMSE_SIC(m_idx,:), 'r-', EbNoVec, BER_ZF_LE(m_idx,:), 'b-', EbNoVec, BER_MMSE_LE(m_idx,:), 'k-', EbNoVec, BER_EGC(m_idx,:), 'b-', EbNoVec, BER_ZF_DF(m_idx,:), 'r-', EbNoVec, BER_MMSE_DF(m_idx,:), 'k-');
        else
            semilogy(EbNoVec, fitBER_MRC, 'r-', EbNoVec, fitBER_MFB, 'k-', EbNoVec, fitBER_ZF_SIC, 'b-', EbNoVec, fitBER_MMSE_SIC, 'r-', EbNoVec, fitBER_ZF_LE, 'b-', EbNoVec, fitBER_MMSE_LE, 'k-', EbNoVec, fitBER_EGC, 'b-', EbNoVec, fitBER_ZF_DF, 'r-', EbNoVec, fitBER_MMSE_DF, 'k-');
        end
        
    end
end
hold off;

m_rx_antennas = '';
for j=1:length(M)
    m_rx_antennas = strcat(m_rx_antennas, sprintf('_%d',M(j)));
end

m_tx_antennas = '';
for j=1:length(K)
    m_tx_antennas = strcat(m_tx_antennas, sprintf('_%d',K(j)));
end

% Get timestamp for saving files.
timeStamp = datestr(now,30);

% Save workspace to MAT-file.
fileName = sprintf('Massive_MU_MIMO_M%s_K%s_flat_fading_various_detectors_%s.mat',m_rx_antennas,m_tx_antennas,timeStamp);
save(fileName);

% Save figure to FIG-file.
fileName = sprintf('Massive_MU_MIMO_M%s_K%s_flat_fading_various_detectors_%s.fig',m_rx_antennas,m_tx_antennas,timeStamp);
savefig(h,fileName);