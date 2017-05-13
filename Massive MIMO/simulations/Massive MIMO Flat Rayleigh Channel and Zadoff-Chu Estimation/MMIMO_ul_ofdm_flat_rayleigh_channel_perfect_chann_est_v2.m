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
EbNoVec = -20:1:-4;                                        % Eb/No in dB.
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
[BER_MRC, BER_ZF_LE, BER_MMSE_LE, BER_MFB] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
h = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink');

% Generate preambles for all the K users. Each PRACH preamble is equivalent to one whole subframe, i.e., 1 ms.
[preambles] = generatePRACHPreamblev2(K);

H = complex(zeros(M,K),zeros(M,K));
%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for idx = 1:length(snr)
    
    linearSnr = 10^(0.1*snr(idx));
    
    nErrs_mrc = 0;
    nErrs_zf_le = 0;
    nErrs_mmse_le = 0;
    nErrs_mfb = 0;
    nBits = 0;
    nBits_mfb = 0;
    ofdm_symbol_number = 0;
    iter = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    while(((nErrs_mrc < nErrors) || (nErrs_zf_le < nErrors) || (nErrs_mmse_le < nErrors) || (nErrs_mfb < nErrors)) && (nBits_mfb < nTotalOfBits))
        
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
            
            %% *********** 1) MRC or Matrix Matched Filter (MMF) receiver *****************
            E_mrc = zeros(modOrd, numSym);
            
            % Apply MRC to received signal.
            G = H';
            r_eq_mrc = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mrc = demodulate(hDemod, r_eq_mrc(jj,:));
                E_mrc(:,jj:K:end) = reshape(demoded_mrc, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ************************ 2) ZF-LE receiver *********************************
            E_zf_le = zeros(modOrd, numSym);
            
            % Apply ZF-LE to received signal.
            H_est = H;
            G = ((H_est'*H_est)^(-1))*H_est';
            r_eq_zf_le = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_zf_le = demodulate(hDemod, r_eq_zf_le(jj,:));
                E_zf_le(:,jj:K:end) = reshape(demoded_zf_le, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ****************************** 3) MMSE-LE receiver *************************
            E_mmse_le = zeros(modOrd, numSym);
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H_est = H;
            G = (H_est'*H_est + (1/linearSnr)*eye(K))^-1;
            G = G*H_est';
            r_eq_mmse = G*ofdm_rx;
            
            % Iterate over all Tx antennas.
            for jj=1:1:K
                demoded_mmse_le = demodulate(hDemod, r_eq_mmse(jj,:));
                E_mmse_le(:,jj:K:end) = reshape(demoded_mmse_le, modOrd, numSym/K);
            end
            % *****************************************************************************
            
            %% ********* 4) Matched Filter Bound (MFB) or Single User detection **********
            
            % Apply MRC to first finger path received signal.
            G = H(:,1)';
            r_eq_mfb = G*ofdm_rx_mfb;
            
            % Demodulate signal.
            E_mfb = demodulate(hDemod, r_eq_mfb);
            % *****************************************************************************
            
            %% -------------------------------- Collect errors ----------------------------
            nErrs_mrc = nErrs_mrc + biterr(msg, E_mrc);
            nErrs_zf_le = nErrs_zf_le + biterr(msg, E_zf_le);
            nErrs_mmse_le = nErrs_mmse_le + biterr(msg, E_mmse_le);
            nErrs_mfb = nErrs_mfb + biterr(msg_mfb, E_mfb);
            
            nBits = nBits + length(msg(:));
            nBits_mfb = nBits_mfb + length(msg_mfb(:));
            fprintf(1,'BER MRC: %f - nErrs_mrc: %d - nBits: %d - iter: %d\n',(nErrs_mrc./nBits),nErrs_mrc,nBits,iter);
            fprintf(1,'BER ZF-LE: %f - nErrs_zf_le: %d - nBits: %d - iter: %d\n',(nErrs_zf_le./nBits),nErrs_zf_le,nBits,iter);
            fprintf(1,'BER MMSE-LE: %f - nErrs_mmse_le: %d - nBits: %d - iter: %d\n',(nErrs_mmse_le./nBits),nErrs_mmse_le,nBits,iter);
            fprintf(1,'BER MFB: %f - nErrs_mfb: %d - nBits_mfb: %d - iter: %d\n',(nErrs_mfb./nBits_mfb),nErrs_mfb,nBits_mfb,iter);
            fprintf(1,'\n');
        else
            ofdm_symbol_number = numSymbInSubframe;
        end
        
        %% --------- Change multipath channel according to its sampling rate. ---------
        if(ofdm_symbol_number==2*numSymbInSubframe)
            ofdm_symbol_number = 0;
        end
        
    end
    
    % Calculate BER for current point
    BER_MRC(idx) = nErrs_mrc./nBits;
    BER_ZF_LE(idx) = nErrs_zf_le./nBits;
    BER_MMSE_LE(idx) = nErrs_mmse_le./nBits;
    BER_MFB(idx) = nErrs_mfb./nBits_mfb;
    
    % Plot results
    semilogy(EbNoVec(1:idx), BER_MRC(1:idx), 'ks', ...
        EbNoVec(1:idx), BER_ZF_LE(1:idx), 'b*', ...
        EbNoVec(1:idx), BER_MMSE_LE(1:idx), 'ro', ...
        EbNoVec(1:idx), BER_MFB(1:idx), 'k');
    legend('MRC', 'ZF-LE', 'MMSE-LE', 'MFB');
    drawnow;
    
end

% Draw the lines.
semilogy(EbNoVec, BER_MRC, 'k-', EbNoVec, BER_ZF_LE, 'b-', EbNoVec, BER_MMSE_LE, 'r-', EbNoVec, BER_MFB, 'k-');
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
    fileName = sprintf('MU_MMIMO_M%s_K%s_flatfading_perfect_chann_est_v%d_%s.fig',m_rx_antennas,m_tx_antennas,script_version,timeStamp);
    savefig(h,fileName);
    
    % Save workspace to MAT-file.
    clear h
    fileName = sprintf('MU_MMIMO_M%s_K%s_flatfading_perfect_chann_est_v%d_%s.mat',m_rx_antennas,m_tx_antennas,script_version,timeStamp);
    save(fileName);
    
end

% Simulation: Massive Multi User MIMO flat-fading and Perfect Channel Estimation
