clear all;close all;clc

script_version = 2;

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                   % Disable/enable figure plotting.

L = 1;                                                              % Number of cells.
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 50;                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
modOrd = 2;                                                         % Constellation size = 2^modOrd.
numSym = K*NFFT;                                                    % Number of symbols, i.e., number of terminals.
NCP = 512;                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = true;                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 0;                                        % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

nTotalOfBits = 1e9;
nErrors = 1e7;

%% ----------------------------------- Set up the simulation -----------------------------------
% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

% Get all bit combinations for ML receiver
bits = de2bi(0:2^(modOrd*K)-1, 'left-msb')';
% Split them per Transmit antenna
b = zeros(K, modOrd, length(bits));
for i = 1:length(bits)
    b(:, :, i) = reshape(bits(:,i), modOrd, K)';
end

% Preallocate variables for speed.
dist = zeros(length(bits), 1);
[BER_MRC_LS, BER_ZF_LS, BER_MMSE_LS, BER_MFB_LS, BER_MRC_MMSE, BER_ZF_MMSE, BER_MMSE_MMSE, BER_MFB_MMSE] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
h = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Uplink MU-MMIMO Flat Rayleigh - M: %d - K: %d - v%d',M,K,script_version);
title(strTitle);

% Generate preambles for all the K users. Each PRACH preamble is equivalent to one whole subframe, i.e., 1 ms.
[preambles] = generatePRACHPreamblev2(K);

H = complex(zeros(M,K),zeros(M,K));
%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for idx = 1:length(snr)
    
    linearSnr = 10^(0.1*snr(idx));
    
    nErrs_mrc_ls = 0;
    nErrs_zf_ls = 0;
    nErrs_mmse_ls = 0;
    nErrs_mfb_ls = 0;
    nErrs_mrc_mmse = 0;
    nErrs_zf_mmse = 0;
    nErrs_mmse_mmse = 0;
    nErrs_mfb_mmse = 0;
    nBits = 0;
    nBits_mfb = 0;
    ofdm_symbol_number = 0;
    iter = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    while(((nErrs_mrc_ls < nErrors) || (nErrs_zf_ls < nErrors) || (nErrs_mmse_ls < nErrors) || (nErrs_mfb_ls < nErrors) || (nErrs_mrc_mmse < nErrors) || (nErrs_zf_mmse < nErrors) || (nErrs_mmse_mmse < nErrors) || (nErrs_mfb_mmse < nErrors)) && (nBits_mfb < nTotalOfBits))
        
        ofdm_symbol_number = ofdm_symbol_number + 1;
        
        if(ofdm_symbol_number > numSymbInSubframe)
            
            %---------- Transmission (UE) ----------
            % Create array of bits to modulate.
            msg = randi(hStr, [0 1], modOrd, numSym);
            
            msg_mfb = msg(:,1:K:end);
            
            % Modulate data.
            source = modulate(hMod, msg);
            
            % Split source among K terminals.
            Tx = reshape(source, K, numel(source)/K); clear source;
            
            % Create OFDM symbol.
            sequence = sqrt(NFFT)*ifft(Tx,NFFT,2);
            
            % Add CP.
            ofdm = [sequence(:,NFFT-NCP+1:end), sequence];
            
            % Make sure the OFDM symbol has unit variance.
            std_dev_vector = (sqrt(diag((ofdm*ofdm')/(NFFT+NCP))));
            std_dev_matrix = diag(1./std_dev_vector);
            signal = std_dev_matrix*ofdm;
        else
            signal = preambles;
        end
        
        %% ---------- Flat Rayleigh Channel plus Noise ----------
        % Flat Rayleigh Fading - independent links. This channel only
        % changes after 24 OFDM symbols, i.e., 2 ms.
        if(ofdm_symbol_number==1)
            H = (1/sqrt(2))*(randn(hStr, M, K) +  1i*randn(hStr, M, K));
        end
        
        % Multi-user.
        x = H*signal;
        
        % Single-user.
        x_mfb = H(:,1)*signal(1,:);
        
        if(add_awgn)
            % Add channel noise power to faded data.
            r = awgn(x, snr(idx), 0, hStr);
            
            % Add channel noise power to faded single user data.
            r_mfb = awgn(x_mfb, snr(idx), 0, hStr);
        else
            r = x;
            r_mfb = x_mfb;
        end
        
        %---------- Reception (base station) ----------
        if(ofdm_symbol_number > numSymbInSubframe)
            
            iter = iter + 1;
            
            % -------- Front End of the receiver --------------------
            % Remove CP.
            rx = r(:,NCP+1:end);
            
            % Retrieve modulation symbols.
            ofdm_rx = (1/sqrt(NFFT))*fft(rx,NFFT,2);
            
            % ----------------- Front End for MFB plotting --------------------
            % Remove CP.
            rx_mfb = r_mfb(:,NCP+1:end);
            
            % Retrieve modulation symbols from first finger path.
            ofdm_rx_mfb = (1/sqrt(NFFT))*fft(rx_mfb,NFFT,2);
            
            % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            %                                  LS Estimation
            % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            
            %% *********** 1) MRC or Matrix Matched Filter (MMF) receiver *****************
            E_mrc_ls = zeros(modOrd, numSym);
            
            % Apply MRC to received signal.
            G = H_estimated_ls';
            r_eq_mrc_ls = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mrc_ls = demodulate(hDemod, r_eq_mrc_ls(jj,:));
                E_mrc_ls(:,jj:K:end) = reshape(demoded_mrc_ls, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ************************ 2) ZF-LE receiver *********************************
            E_zf_ls = zeros(modOrd, numSym);
            
            % Apply ZF-LE to received signal.
            H_est = H_estimated_ls;
            G = ((H_est'*H_est)^(-1))*H_est';
            r_eq_zf_le_ls = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_zf_le_ls = demodulate(hDemod, r_eq_zf_le_ls(jj,:));
                E_zf_ls(:,jj:K:end) = reshape(demoded_zf_le_ls, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ****************************** 3) MMSE-LE receiver *************************
            E_mmse_ls = zeros(modOrd, numSym);
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H_est = H_estimated_ls;
            G = (H_est'*H_est + (1/linearSnr)*eye(K))^-1;
            G = G*H_est';
            r_eq_mmse_ls = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mmse_ls = demodulate(hDemod, r_eq_mmse_ls(jj,:));
                E_mmse_ls(:,jj:K:end) = reshape(demoded_mmse_ls, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ********* 4) Matched Filter Bound (MFB) or Single User detection **********
            
            % Apply MRC to first finger path received signal.
            G = H_estimated_ls(:,1)';
            r_eq_mfb_ls = G*ofdm_rx_mfb;
            
            % Demodulate signal.
            E_mfb_ls = demodulate(hDemod, r_eq_mfb_ls);
            % *****************************************************************************
            
            % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            %                                  MMSE Estimation
            % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            
            %% *********** 1) MRC or Matrix Matched Filter (MMF) receiver *****************
            E_mrc_mmse = zeros(modOrd, numSym);
            
            % Apply MRC to received signal.
            G = H_estimated_mmse';
            r_eq_mrc_mmse = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mrc_mmse = demodulate(hDemod, r_eq_mrc_mmse(jj,:));
                E_mrc_mmse(:,jj:K:end) = reshape(demoded_mrc_mmse, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ************************ 2) ZF-LE receiver *********************************
            E_zf_mmse = zeros(modOrd, numSym);
            
            % Apply ZF-LE to received signal.
            H_est = H_estimated_mmse;
            G = ((H_est'*H_est)^(-1))*H_est';
            r_eq_zf_mmse = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_zf_mmse = demodulate(hDemod, r_eq_zf_mmse(jj,:));
                E_zf_mmse(:,jj:K:end) = reshape(demoded_zf_mmse, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ****************************** 3) MMSE-LE receiver *************************
            E_mmse_mmse = zeros(modOrd, numSym);
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H_est = H_estimated_mmse;
            G = (H_est'*H_est + (1/linearSnr)*eye(K))^-1;
            G = G*H_est';
            r_eq_mmse_mmse = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mmse_mmse = demodulate(hDemod, r_eq_mmse_mmse(jj,:));
                E_mmse_mmse(:,jj:K:end) = reshape(demoded_mmse_mmse, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ********* 4) Matched Filter Bound (MFB) or Single User detection **********
            
            % Apply MRC to first finger path received signal.
            G = H_estimated_mmse(:,1)';
            r_eq_mfb_mmse = G*ofdm_rx_mfb;
            
            % Demodulate signal.
            E_mfb_mmse = demodulate(hDemod, r_eq_mfb_mmse);
            % *****************************************************************************
            
            %% -------------------------------- Collect errors ----------------------------
            nErrs_mrc_ls = nErrs_mrc_ls + biterr(msg, E_mrc_ls);
            nErrs_zf_ls = nErrs_zf_ls + biterr(msg, E_zf_ls);
            nErrs_mmse_ls = nErrs_mmse_ls + biterr(msg, E_mmse_ls);
            nErrs_mfb_ls = nErrs_mfb_ls + biterr(msg_mfb, E_mfb_ls);
            nErrs_mrc_mmse = nErrs_mrc_mmse + biterr(msg, E_mrc_mmse);
            nErrs_zf_mmse = nErrs_zf_mmse + biterr(msg, E_zf_mmse);
            nErrs_mmse_mmse = nErrs_mmse_mmse + biterr(msg, E_mmse_mmse);
            nErrs_mfb_mmse = nErrs_mfb_mmse + biterr(msg_mfb, E_mfb_mmse);
            
            nBits = nBits + length(msg(:));
            nBits_mfb = nBits_mfb + length(msg_mfb(:));
            fprintf(1,'BER LS/MRC: %f - nErrs_mrc_ls: %d - nBits: %d - iter: %d\n',(nErrs_mrc_ls./nBits),nErrs_mrc_ls,nBits,iter);
            fprintf(1,'BER MMSE/MRC: %f - nErrs_mrc_mmse: %d - nBits: %d - iter: %d\n',(nErrs_mrc_mmse./nBits),nErrs_mrc_mmse,nBits,iter);
            fprintf(1,'BER LS/ZF: %f - nErrs_zf_ls: %d - nBits: %d - iter: %d\n',(nErrs_zf_ls./nBits),nErrs_zf_ls,nBits,iter);
            fprintf(1,'BER MMSE/ZF: %f - nErrs_zf_mmse: %d - nBits: %d - iter: %d\n',(nErrs_zf_mmse./nBits),nErrs_zf_mmse,nBits,iter);
            fprintf(1,'BER LS/MMSE: %f - nErrs_mmse_ls: %d - nBits: %d - iter: %d\n',(nErrs_mmse_ls./nBits),nErrs_mmse_ls,nBits,iter);
            fprintf(1,'BER MMSE/MMSE: %f - nErrs_mmse_mmse: %d - nBits: %d - iter: %d\n',(nErrs_mmse_mmse./nBits),nErrs_mmse_mmse,nBits,iter);
            fprintf(1,'BER LS/MFB: %f - nErrs_mfb_ls: %d - nBits_mfb: %d - iter: %d\n',(nErrs_mfb_ls./nBits_mfb),nErrs_mfb_ls,nBits_mfb,iter);
            fprintf(1,'BER MMSE/MFB: %f - nErrs_mfb_mmse: %d - nBits_mfb: %d - iter: %d\n',(nErrs_mfb_mmse./nBits_mfb),nErrs_mfb_mmse,nBits_mfb,iter);
            fprintf(1,'\n');
        else
            [ID, TA, H_estimated_ls] = detectPreambleIDAndTAv5(r, M, K);
            
            var_channel = 1;
            var_noise = 1/linearSnr;
            [ID, TA, H_estimated_mmse] = detectPreambleIDAndTAv6(r, M, K, var_channel, var_noise);
            
            ofdm_symbol_number = numSymbInSubframe;
            
%             error_ls = sum(sum(abs(H-H_estimated_ls).^2))/(M*K);
%             error_mmse = sum(sum(abs(H-H_estimated_mmse).^2))/(M*K);
%             
%             fprintf(1,'error_ls: %d\n',error_ls);
%             fprintf(1,'error_mmse: %d\n',error_mmse);
%             fprintf(1,'\n');
        end
        
        %% --------- Change multipath channel according to its sampling rate. ---------
        if(ofdm_symbol_number==2*numSymbInSubframe)
            ofdm_symbol_number = 0;
        end
        
    end
    
    % Calculate BER for current point
    BER_MRC_LS(idx) = nErrs_mrc_ls./nBits;
    BER_ZF_LS(idx) = nErrs_zf_ls./nBits;
    BER_MMSE_LS(idx) = nErrs_mmse_ls./nBits;
    BER_MFB_LS(idx) = nErrs_mfb_ls./nBits_mfb;
    BER_MRC_MMSE(idx) = nErrs_mrc_mmse./nBits;
    BER_ZF_MMSE(idx) = nErrs_zf_mmse./nBits;
    BER_MMSE_MMSE(idx) = nErrs_mmse_mmse./nBits;
    BER_MFB_MMSE(idx) = nErrs_mfb_mmse./nBits_mfb;
    
    % Plot results
    semilogy(EbNoVec(1:idx), BER_MRC_LS(1:idx), 'ks', ...
        EbNoVec(1:idx), BER_MRC_MMSE(1:idx), 'ko', ...
        EbNoVec(1:idx), BER_ZF_LS(1:idx), 'bs', ...
        EbNoVec(1:idx), BER_ZF_MMSE(1:idx), 'bo', ...
        EbNoVec(1:idx), BER_MMSE_LS(1:idx), 'rs', ...
        EbNoVec(1:idx), BER_MMSE_MMSE(1:idx), 'ro', ...
        EbNoVec(1:idx), BER_MFB_LS(1:idx), 'gs', ...
        EbNoVec(1:idx), BER_MFB_MMSE(1:idx), 'go');
    legend('LS/MRC', 'MMSE/MRC', 'LS/ZF', 'MMSE/ZF', 'LS/MMSE', 'MMSE/MMSE', 'LS/MFB', 'MMSE/MFB');
    drawnow;
    
end

% Draw the lines.
semilogy(EbNoVec, BER_MRC_LS, 'k-', ...
    EbNoVec, BER_MRC_MMSE, 'k-', ...
    EbNoVec, BER_ZF_LS, 'b-', ...
    EbNoVec, BER_ZF_MMSE, 'b-', ...
    EbNoVec, BER_MMSE_LS, 'r-', ...
    EbNoVec, BER_MMSE_MMSE, 'r-', ...
    EbNoVec, BER_MFB_LS, 'g-', ...
    EbNoVec, BER_MFB_MMSE, 'g-');
hold off;

if(add_awgn)
    
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
    
    % Save figure to FIG-file.
    fileName = sprintf('M_MU_MIMO_M%s_K%s_flatfading_zc_est_chann_v%d_%s_.fig',m_rx_antennas,m_tx_antennas,script_version,timeStamp);
    savefig(h,fileName);
    
    % Save workspace to MAT-file.
    clear h
    fileName = sprintf('M_MU_MIMO_M%s_K%s_flatfading_zc_est_chann_v%d_%s.mat',m_rx_antennas,m_tx_antennas,script_version,timeStamp);
    save(fileName);
    
end

% Simulation: Massive Multi User MIMO flat-fading and Zadoff-Chu Estimated Channel
