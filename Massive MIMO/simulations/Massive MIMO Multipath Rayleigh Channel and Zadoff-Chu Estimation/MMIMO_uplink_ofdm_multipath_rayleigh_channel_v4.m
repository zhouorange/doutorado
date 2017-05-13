clear all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                   % Disable/enable figure plotting.

L = 1;                                                              % Number of cells.
K = 5;                                                              % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 50;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
modOrd = 2;                                                         % Constellation size = 2^modOrd.
numSym = K*NFFT;                                                    % Number of symbols, i.e., number of terminals.
NCP = 512;                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
Ts = 1/(15000*NFFT);                                                % System Sampleing Rate.
numSymbInSubframe = 12;                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

EbNoVec = -20:1:0;                                                  % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e7;
nErrors = 100000;
debug = false;

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 0.4]*1e-6;                                               % Delay in microseconds.
gain  = [0 0];                                                      % Gain in dB (indoor).
numPaths = length(delay);                                           % Number of paths per channel.
totalNumPaths = M*K*numPaths;                                       % Total number of paths between the various antennas and Base Stations.

%% ----------------------------------- Setup Channel -----------------------------------
% Vetor de ganhos.
pos = round(delay/Ts)+1;                                            % Effective position of taps within the vector.
g = zeros(1, round(delay(end)/Ts)+1);                               % +1 is used to include the delay 0 into the vector.
for n = 1:length(delay)
    g( pos(n) ) = 10^( gain(n)/10 );
end

fc = 1e9;                                                           % Carrier Freq. in MHz.
c = 3e8;                                                            % Light speed in m/s.
v = 30;                                                             % Speed in m/s.
Fs = 30.72e6;                                                       % Sampling Freq. in MHz.
Ts = 1/Fs;                                                          % Sampling period in seconds.

fd = (v*fc)/c;                                                      % Doppler frequency in Hz.

Pd = 0;                                                             % Relative power in dB.

% Parameters used for generating the Mulipath Masive MIMO Channel.
Fs_chann = 500;                                                     % Channel Sampling Rate in Hz. The channel is sampled every 2 ms. (Periodo de amostragem do canal ~1 amostra / Slot OFDM)
Ts_chann = 1/Fs_chann;
N_chann = 256;                                                      % Number of samples used to sample the channel. Duration of the channel in  number of samples.
delta_f_chann = Fs_chann/N_chann;                                   % in Hz.
f = -Fs_chann/2:delta_f_chann:Fs_chann/2;
idx = find(f<fd);
f = f(idx);
idx = find(f>-fd);
f = f(idx);
f = f.';

LS = length(f);
S = 1/pi/fd./sqrt(1 - (f/fd).^2) * 10^(Pd/10);
S = S * LS / sum(S);                                                % Energy nomalization.
S1 = S;
S = [S((LS+1)/2:LS); zeros(N_chann-LS,1); S(1:(LS-1)/2)];           % Interpolation of the Doppler Spectrum by N_chann/LS.

% ************************* Generate Multipath Massive MIMO Channel. **********************
rng(55);
x = [(randn((LS-1)/2+1,totalNumPaths,'double') + 1i*randn((LS-1)/2+1,totalNumPaths,'double')); zeros(N_chann-LS,totalNumPaths); (randn((LS-1)/2,totalNumPaths,'double')) + 1i*randn((LS-1)/2,totalNumPaths,'double')]/sqrt(2);
ch = ifft(x .* repmat(sqrt(S),1,totalNumPaths)) * N_chann / sqrt(LS);

% Plot doppler spectrum figures.
if(plot_fig)
    figure;
    plot(f, abs(S1));
    
    figure;
    plot(S)
    
    figure;
    Tsd = 1/Fs_chann;
    TNd = Tsd*N_chann; % Channel duration in seconds.
    plot( (0:N_chann-1)*Tsd, 10*log10(abs(ch)) );
    xlabel('Time (s)')
    ylabel('Power (dB)');
    axis([0 TNd -30 10]);
end

% **************************** Cell Layout. *******************************
% Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
radius = cellHole + (cellRadius-cellHole).*rand(1,K);

% Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
angle = 2*pi*rand(1,K);

% Plot the position of each terminal inside the cell.
figure;
plot(0, 0, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]);
hold on

Circlex = cellRadius*cos(2*pi*(0:99)/100);
Circley = cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

Circlex = cellHole*cos(2*pi*(0:99)/100);
Circley = cellHole*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'r');

UEx = radius.*cos(angle);
UEy = radius.*sin(angle);
plot(UEx, UEy, 'b*');

grid on;
hold off;

% Calculate path-loss for all users, in meters.
path_loss = radius.^gamma;

% ************ Apply power delay profile and large-scale fading to the multiptah channel matrix, i.e., small-scale fading. ************
for idx_ch = 1 : N_chann
    
    if(power_ctrl_enabled)
        largeScaleFading = 1;
    else
        % Calculate shadowing for each one of different channels.
        shadowing_att = lognrnd(0,sshadow,1,K);
        % Large scale fading calculated according to Marzetta.
        largeScaleFading = repmat(sqrt(shadowing_att./path_loss), M, 1);
    end
    
    % Atualizacao das matrizes de canal para a celula alvo (l = 1): existe uma matriz de canal por percurso do terminal.
    G = reshape(ch(idx_ch,:), M, K, length(pos));
    G(:,:,1) = (g(pos(1)) * G(:,:,1)) .* largeScaleFading;
    for k = 2:length(pos)
        G(:,:,k) = (g(pos(k)) * G(:,:,k)) .* largeScaleFading;
    end
    ch(idx_ch,:) = reshape(G, 1, totalNumPaths);
end

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
[BER_ZF_SIC, BER_MMSE_SIC, BER_ML, BER_MRC, BER_ZF_LE, BER_MMSE_LE, BER_EGC, BER_ZF_DF, BER_MMSE_DF, BER_MFB] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
h = gcf; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink');

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for idx = 1:length(snr)
    
    linearSnr = 10^(0.1*snr(idx));
    
    nErrs_zf_sic = 0;
    nErrs_mmse_sic = 0;
    nErrs_ml = 0;
    nErrs_mrc = 0;
    nErrs_zf_le = 0;
    nErrs_mmse_le = 0;
    nErrs_egc = 0;
    nErrs_zf_df = 0;
    nErrs_mmse_df = 0;
    nErrs_mfb = 0;
    nBits = 0;
    nBits_mfb = 0;
    ofdm_symbol_number = 0;
    idx_ch = 1;
    iter = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    while(((nErrs_zf_sic < nErrors) || (nErrs_mmse_sic < nErrors) || (nErrs_ml < nErrors) || (nErrs_mrc < nErrors) || (nErrs_zf_le < nErrors) || (nErrs_mmse_le < nErrors) || (nErrs_egc < nErrors) || (nErrs_zf_df < nErrors) || (nErrs_mmse_df < nErrors) || (nErrs_mfb < nErrors)) && (nBits_mfb < nTotalOfBits))
        
        ofdm_symbol_number = ofdm_symbol_number + 1;
        
        iter = iter + 1;
        
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
        ofdm = std_dev_matrix*ofdm;
        
        %---------- Multipath Channel plus Noise ----------
        H = reshape(ch(idx_ch,:), M, K, length(pos));
        
        x = H(:,:,1)*ofdm;
        for k = 2:length(pos)
            aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*ofdm(:,1:end-(pos(k)-1))];
            x = x + aux;
        end
        
        % Single User, used for plotting the Matched Matrix Bound (MFB).
        x_mfb = [H(:,1,1)*ofdm(1,:) complex(zeros(M,(pos(length(pos))-1)),zeros(M,(pos(length(pos))-1)))];
        for k = 2:length(pos)
            aux_mfb = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,1,k)*ofdm(1,:)];
            x_mfb = x_mfb + aux_mfb;
        end
        
        % Add channel noise power to faded data.
        r = awgn(x, snr(idx), 0, hStr);
        
        % Add channel noise power to faded single user data.
        r_mfb = awgn(x_mfb, snr(idx), 0, hStr);
        
        %---------- Reception (base station) ----------
        
        % Remove CP.
        rx = r(:,NCP+1:end);
        
        % Retrieve modulation symbols.
        ofdm_rx = (1/sqrt(NFFT))*fft(rx,NFFT,2);
        
        %% *********** 1) Zero-Forcing with Optimally Ordered SIC receiver ************
        E_zf_sic = zeros(modOrd, numSym); k = zeros(K, 1);
        
        % Move OFDM symbol to aux variable.
        ofdm_rx_zf_sic = ofdm_rx;
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Initialization
        G = pinv(H_estimated);
        [val, k0] = min(sum(abs(G).^2,2));
        
        % Start Zero-Forcing Nulling Loop.
        for n = 1:K
            % Find best transmitter signal using minimum norm.
            k(n) = k0;
            
            % Select Weight vector for best transmitter signal.
            w = G(k(n),:);
            
            % Calculate output for transmitter n.
            y = w * ofdm_rx_zf_sic;
            
            % Demodulate bitstream.
            demoded_zf_sic = demodulate(hDemod, y);
            E_zf_sic(:, k(n):K:end) = reshape(demoded_zf_sic, modOrd, numSym/K);
            
            % Subtract effect of the transmitter n from received signal.
            z = modulate(hMod, demodulate(hDemod, y));
            ofdm_rx_zf_sic = ofdm_rx_zf_sic - H_estimated(:, k(n))*z;
            
            % Adjust channel estimate matrix for next minimum norm search.
            H_estimated(:, k(n)) = zeros(M, 1);
            G = pinv(H_estimated);
            for aa = 1:n
                G(k(aa), :) = inf;
            end
            [val, k0] = min(sum(abs(G).^2,2));
        end
        % *****************************************************************************
        
        %% ************** 2) MMSE with Optimally Ordered SIC receiver *****************
        E_mmse_sic = zeros(modOrd, numSym); k = zeros(K, 1);
        
        % Move OFDM symbol to aux variable.
        ofdm_rx_mmse_sic = ofdm_rx;
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Initialization
        G = (H_estimated*H_estimated' + (1/linearSnr)*eye(M))^-1;
        G = H_estimated'*G;
        [val, k0] = min(sum(abs(G).^2,2));
        % Start MMSE Nulling Loop
        for n = 1:K
            % Find best transmitter signal using Min Norm
            k(n) = k0;
            
            % Select Weight vector for best transmitter signal
            w = G(k(n),:);
            
            % Calculate output for transmitter n and demodulate bitstream
            y = w * ofdm_rx_mmse_sic;
            E_mmse_sic(:, k(n):K:end) = reshape(demodulate(hDemod, y), modOrd, numSym/K);
            
            % Subtract effect of the transmitter n from received signal
            z = modulate(hMod, demodulate(hDemod, y));
            ofdm_rx_mmse_sic = ofdm_rx_mmse_sic - H_estimated(:, k(n))*z;
            
            % Adjust channel estimate matrix for next min Norm search
            H_estimated(:, k(n)) = zeros(M, 1);
            G = (H_estimated*H_estimated' + (1/linearSnr)*eye(M))^-1;
            G = H_estimated'*G;
            for aa = 1:n
                G(k(aa), :) = inf;
            end
            [val, k0] = min(sum(abs(G).^2,2));
        end
        % *****************************************************************************
        
        %% ************** 3) Joint Maximum Likelihood (JML) receiver ******************
        E_ml = zeros(modOrd, numSym);
        
        % Move OFDM symbol to aux variable.
        ofdm_rx_ml = ofdm_rx;
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        for j=1:NFFT
            for i = 1:2^(modOrd*K)
                % Signal constellation for each bit combination.
                sig = modulate(hMod, b(:, :, i)').';
                
                % Distance metric for each constellation.
                dist(i) = sum(abs(ofdm_rx_ml(:,j) - H_estimated*sig).^2);
            end
            % Get the minimum.
            [notUsed, val] = min(dist);
            % detected bits.
            E_ml(:,((j-1)*K)+1:j*K) = b(:,:,val)';
        end
        % *****************************************************************************
        
        %% *********** 4) MRC or Matrix Matched Filter (MMF) receiver *****************
        E_mrc = zeros(modOrd, numSym);
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Apply MRC to received signal.
        G = H_estimated';
        r_mrc = G*ofdm_rx;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_mrc = demodulate(hDemod, r_mrc(jj,:));
            E_mrc(:,jj:K:end) = reshape(demoded_mrc, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ************************ 5) ZF-LE receiver *********************************
        E_zf_le = zeros(modOrd, numSym);
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Apply ZF-LE to received signal.
        G = ((H_estimated'*H_estimated)^(-1))*H_estimated';
        r_zf_le = G*ofdm_rx;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_zf_le = demodulate(hDemod, r_zf_le(jj,:));
            E_zf_le(:,jj:K:end) = reshape(demoded_zf_le, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ****************************** 6) MMSE-LE receiver *************************
        E_mmse_le = zeros(modOrd, numSym);
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Apply MMSE-LE to received signal.
        G = (H_estimated*H_estimated' + (1/linearSnr)*eye(M))^-1;
        G = H_estimated'*G;
        r_mmse_le = G*ofdm_rx;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_mmse_le = demodulate(hDemod, r_mmse_le(jj,:));
            E_mmse_le(:,jj:K:end) = reshape(demoded_mmse_le, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% *********************** 8) Equal Gain Combining (EGC) **********************
        E_egc = zeros(modOrd, numSym);
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Apply EGC to received signal.
        G = exp(-1*1i*angle(H_estimated)).';
        
        % Remove the phase of the channel.
        r_egc = G*ofdm_rx;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_egc = demodulate(hDemod, r_egc(jj,:));
            E_egc(:,jj:K:end) = reshape(demoded_egc, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ************* 9) Zero-Forcing Decision Feedback receiver (ZF-DF) ************
        E_zf_df = zeros(modOrd, numSym);
        
        % Move OFDM symbol to aux variable.
        ofdm_rx_zf_df = ofdm_rx;
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Initialization
        G = pinv(H_estimated);
        % Start Zero-Forcing Nulling Loop.
        for n = 1:K
            % Select Weight vector for best transmitter signal.
            w = G(n,:);
            
            % Calculate output for transmitter n.
            y = w * ofdm_rx_zf_df;
            
            % Demodulate bitstream.
            demoded_zf_df = demodulate(hDemod, y);
            E_zf_df(:, n:K:end) = reshape(demoded_zf_df, modOrd, numSym/K);
            
            % Subtract effect of the transmitter n from received signal.
            z = modulate(hMod, demodulate(hDemod, y));
            ofdm_rx_zf_df = ofdm_rx_zf_df - H_estimated(:, n)*z;
        end
        % *****************************************************************************
        
        %% ************** 10) MMSE-Decision Feedback receiver (MMSE-DF) ***************
        E_mmse_df = zeros(modOrd, numSym);
        
        % Move OFDM symbol to aux variable.
        ofdm_rx_mmse_df = ofdm_rx;
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Initialization.
        G = (H_estimated*H_estimated' + (1/linearSnr)*eye(M))^-1;
        G = H_estimated'*G;
        % Start MMSE Nulling Loop.
        for n = 1:K
            % Select Weight vector for best transmitter signal.
            w = G(n,:);
            
            % Calculate output for transmitter n and demodulate bitstream.
            y = w * ofdm_rx_mmse_df;
            E_mmse_df(:, n:K:end) = reshape(demodulate(hDemod, y), modOrd, numSym/K);
            
            % Subtract effect of the transmitter n from received signal.
            z = modulate(hMod, demodulate(hDemod, y));
            ofdm_rx_mmse_df = ofdm_rx_mmse_df - H_estimated(:, n)*z;
        end
        % *****************************************************************************
        
        %% ********* 11) Matched Filter Bound (MFB) or Single User detection **********
        
        % Assume perfect channel estimation.
        H_estimated = H(:,:,1);
        
        % Remove CP.
        rx_mfb = r_mfb(:,NCP+1:end);
        
        % Retrieve modulation symbols.
        ofdm_rx_mfb = (1/sqrt(NFFT))*fft(rx_mfb,NFFT,2);
        
        % Apply MRC to received signal.
        G = H_estimated(:,1)';
        r_eq_mgb = G*ofdm_rx_mfb;
        
        % Demodulate signal.
        E_mfb = demodulate(hDemod, r_eq_mgb);
        % *****************************************************************************
        
        %% --------- Change multipath channel according to its sampling rate. ---------
        if(ofdm_symbol_number==2*numSymbInSubframe)
            ofdm_symbol_number = 0;
            idx_ch = idx_ch + 1;
            if(idx_ch > N_chann)
                idx_ch = 1;
            end
        end
        
        %% -------------------------------- Collect errors ----------------------------
        nErrs_zf_sic = nErrs_zf_sic + biterr(msg, E_zf_sic);
        nErrs_mmse_sic = nErrs_mmse_sic + biterr(msg, E_mmse_sic);
        nErrs_ml = nErrs_ml + biterr(msg, E_ml);
        nErrs_mrc = nErrs_mrc + biterr(msg, E_mrc);
        nErrs_zf_le = nErrs_zf_le + biterr(msg, E_zf_le);
        nErrs_mmse_le = nErrs_mmse_le + biterr(msg, E_mmse_le);
        nErrs_egc = nErrs_egc + biterr(msg, E_egc);
        nErrs_zf_df = nErrs_zf_df + biterr(msg, E_zf_df);
        nErrs_mmse_df = nErrs_mmse_df + biterr(msg, E_mmse_df);
        nErrs_mfb = nErrs_mfb + biterr(msg_mfb, E_mfb);
        
        nBits = nBits + length(msg(:));
        nBits_mfb = nBits_mfb + length(msg_mfb(:));
        fprintf(1,'BER ZF-SIC: %f - nErrs_zf_sic: %d - nBits: %d - iter: %d\n',(nErrs_zf_sic./nBits),nErrs_zf_sic,nBits,iter);
        fprintf(1,'BER MMSE-SIC: %f - nErrs_mmse_sic: %d - nBits: %d - iter: %d\n',(nErrs_mmse_sic./nBits),nErrs_mmse_sic,nBits,iter);
        fprintf(1,'BER MRC: %f - nErrs_mrc: %d - nBits: %d - iter: %d\n',(nErrs_mrc./nBits),nErrs_mrc,nBits,iter);
        fprintf(1,'BER ZF-LE: %f - nErrs_zf_le: %d - nBits: %d - iter: %d\n',(nErrs_zf_le./nBits),nErrs_zf_le,nBits,iter);
        fprintf(1,'BER MMSE-LE: %f - nErrs_mmse_le: %d - nBits: %d - iter: %d\n',(nErrs_mmse_le./nBits),nErrs_mmse_le,nBits,iter);
        fprintf(1,'BER EGC: %f - nErrs_egc: %d - nBits: %d - iter: %d\n',(nErrs_egc./nBits),nErrs_egc,nBits,iter);
        fprintf(1,'BER ZF-DF: %f - nErrs_zf_df: %d - nBits: %d - iter: %d\n',(nErrs_zf_df./nBits),nErrs_zf_df,nBits,iter);
        fprintf(1,'BER MMSE-DF: %f - nErrs_mmse_df: %d - nBits: %d - iter: %d\n',(nErrs_mmse_df./nBits),nErrs_mmse_df,nBits,iter);
        fprintf(1,'BER ML: %f - nErrs_ml: %d - nBits: %d - iter: %d\n',(nErrs_ml./nBits),nErrs_ml,nBits,iter);
        fprintf(1,'BER MFB: %f - nErrs_mfb: %d - nBits_mfb: %d - iter: %d\n',(nErrs_mfb./nBits_mfb),nErrs_mfb,nBits_mfb,iter);
        fprintf(1,'\n');
        
    end
    
    % Calculate BER for current point
    BER_ZF_SIC(idx) = nErrs_zf_sic./nBits;
    BER_MMSE_SIC(idx) = nErrs_mmse_sic./nBits;
    BER_MRC(idx) = nErrs_mrc./nBits;
    BER_ZF_LE(idx) = nErrs_zf_le./nBits;
    BER_MMSE_LE(idx) = nErrs_mmse_le./nBits;
    BER_EGC(idx) = nErrs_egc./nBits;
    BER_ZF_DF(idx) = nErrs_zf_df./nBits;
    BER_MMSE_DF(idx) = nErrs_mmse_df./nBits;
    BER_ML(idx) = nErrs_ml./nBits;
    BER_MFB(idx) = nErrs_mfb./nBits_mfb;
    
    % Plot results
    semilogy(EbNoVec(1:idx), BER_ZF_SIC(1:idx), 'r*', ...
        EbNoVec(1:idx), BER_MMSE_SIC(1:idx), 'bo', ...
        EbNoVec(1:idx), BER_ML(1:idx), 'gs', ...
        EbNoVec(1:idx), BER_MRC(1:idx), 'ks', ...
        EbNoVec(1:idx), BER_ZF_LE(1:idx), 'b*', ...
        EbNoVec(1:idx), BER_MMSE_LE(1:idx), 'ko', ...
        EbNoVec(1:idx), BER_EGC(1:idx), 'go', ...
        EbNoVec(1:idx), BER_ZF_DF(1:idx), 'rs', ...
        EbNoVec(1:idx), BER_MMSE_DF(1:idx), 'bs', ...
        EbNoVec(1:idx), BER_MFB(1:idx), 'ys');
    legend('ZF-SIC', 'MMSE-SIC', 'ML', 'MRC', 'ZF-LE', 'MMSE-LE', 'EGC', 'ZF-DF', 'MMSE-DF', 'MFB');
    drawnow;
    
end

% Draw the lines
semilogy(EbNoVec, BER_ZF_SIC, 'r-', EbNoVec, BER_MMSE_SIC, 'b-', ...
    EbNoVec, BER_ML, 'g-', EbNoVec, BER_MRC, 'k-', EbNoVec, BER_ZF_LE, 'b-', EbNoVec, BER_MMSE_LE, 'k-', EbNoVec, BER_EGC, 'g-', EbNoVec, BER_ZF_DF, 'r-', EbNoVec, BER_MMSE_DF, 'b-', EbNoVec, BER_MFB, 'y-');
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
fileName = sprintf('Massive_MU_MIMO_M%s_K%s_multipath_fading_various_detectors_%s.mat',m_rx_antennas,m_tx_antennas,timeStamp);
save(fileName);

% Save figure to FIG-file.
fileName = sprintf('Massive_MU_MIMO_M%s_K%s_multipath_fading_various_detectors_%s.fig',m_rx_antennas,m_tx_antennas,timeStamp);
savefig(h,fileName);
