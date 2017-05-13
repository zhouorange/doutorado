clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                   % Disable/enable figure plotting.

L = 1;                                                              % Number of cells.
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
modOrd = 2;                                                         % Constellation size = 2^modOrd.
numSym = K*NFFT;                                                    % Number of symbols, i.e., number of terminals.
NCP = 512;                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = false;                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 1500;%-20:1:-4;                                                 % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e9;
nErrors = 10000000;
debug = false;

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 0.977]*1e-6;                                             % Delay in microseconds.
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
[BER_MRC, BER_ZF_LE, BER_MMSE_LE, BER_MFB] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
%h = gcf; 
h = figure;
grid on; hold on;
set(gca,'yscale','log','xlim',[(EbNoVec(1)-0.01), EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink');

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
    idx_ch = 1;
    iter = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    while(((nErrs_mrc < nErrors) || (nErrs_zf_le < nErrors) || (nErrs_mmse_le < nErrors) || (nErrs_mfb < nErrors)) && (nBits_mfb < nTotalOfBits))
        
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
        
        x = [H(:,:,1)*ofdm complex(zeros(M,(pos(length(pos))-1)),zeros(M,(pos(length(pos))-1)))];
        for k = 2:length(pos)
            aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*ofdm];
            x = x + aux;
        end
        
        % Single User, used for plotting the Matched Matrix Bound (MFB).
        x_mfb = [H(:,1,1)*ofdm(1,:) complex(zeros(M,(pos(length(pos))-1)),zeros(M,(pos(length(pos))-1)))];
        for k = 2:length(pos)
            aux_mfb = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,1,k)*ofdm(1,:)];
            x_mfb = x_mfb + aux_mfb;
        end
        
        if(add_awgn)
            % Add channel noise power to faded data.
            r = awgn(x, snr(idx), 0, hStr);
            
            % Add channel noise power to faded single user data.
            r_mfb = awgn(x_mfb, snr(idx), 0, hStr);
        else
            r = x;
            r_mfb = x_mfb;
        end
        
        %-------------------- Reception (base station) --------------------
        % @@@@@@@@@@ Assume perfect channel estimation. @@@@@@@@@@
        
        % -------- Front End of the Rake-like receiver --------------------   
        % Get fisrt path signal, i.e., the most powerful signal, i.e., the direct path.
        r_first_finger_path = r(:,1:end-(pos(length(pos))-1));
                
        % get the second path signal.
        r_second_finger_path = r(:,pos(length(pos)):end);
        
        % Remove CP from first finger path.
        rx_finger1 = r_first_finger_path(:,NCP+1:end);
        
        % Remove CP from second finger path.
        rx_finger2 = r_second_finger_path(:,NCP+1:end);
        
        % Retrieve modulation symbols from first finger path.
        ofdm_rx_finger1 = (1/sqrt(NFFT))*fft(rx_finger1,NFFT,2);
        
        % Retrieve modulation symbols from second finger path.
        ofdm_rx_finger2 = (1/sqrt(NFFT))*fft(rx_finger2,NFFT,2);
        
        % ----------------- Front End for MFB plotting --------------------
        % Get fisrt path signal, i.e., the most powerful signal, i.e., the direct path.
        r_first_finger_path_mfb = r_mfb(:,1:end-(pos(length(pos))-1));
                
        % get the second path signal.
        r_second_finger_path_mfb = r_mfb(:,pos(length(pos)):end);
        
        % Remove CP from direct path.
        rx_finger1_mfb = r_first_finger_path_mfb(:,NCP+1:end);
        
        % Remove CP from second path.
        rx_finger2_mfb = r_second_finger_path_mfb(:,NCP+1:end);
        
        % Retrieve modulation symbols from first finger path.
        ofdm_rx_finger1_mfb = (1/sqrt(NFFT))*fft(rx_finger1_mfb,NFFT,2);
        
        % Retrieve modulation symbols from second finger path.
        ofdm_rx_finger2_mfb = (1/sqrt(NFFT))*fft(rx_finger2_mfb,NFFT,2);
        
        %% *********** 1) MRC or Matrix Matched Filter (MMF) receiver *****************
        E_mrc = zeros(modOrd, numSym);
        
        % Apply MRC to first finger path received signal.
        G1 = H(:,:,1)';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MRC to second finger path received signal.
        G2 = H(:,:,2)';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mrc = (r_eq1 + r_eq2)/2;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_mrc = demodulate(hDemod, r_eq_sum_mrc(jj,:));
            E_mrc(:,jj:K:end) = reshape(demoded_mrc, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ************************ 2) ZF-LE receiver *********************************
        E_zf_le = zeros(modOrd, numSym);
        
        % Apply ZF-LE to first finger path received signal.
        G1 = ((H(:,:,1)'*H(:,:,1))^(-1))*H(:,:,1)';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply ZF-LE to second path received signal.
        G2 = ((H(:,:,2)'*H(:,:,2))^(-1))*H(:,:,2)';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_zf_le = (r_eq1 + r_eq2)/2;
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_zf_le = demodulate(hDemod, r_eq_sum_zf_le(jj,:));
            E_zf_le(:,jj:K:end) = reshape(demoded_zf_le, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ****************************** 3) MMSE-LE receiver *************************
        E_mmse_le = zeros(modOrd, numSym);
        
        % Apply MMSE-LE to fisrt finger path received signal.
        H1_est = H(:,:,1);
        G1 = (H1_est'*H1_est + (1/linearSnr)*eye(K))^-1;
        G1 = G1*H1_est';   
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MMSE-LE to second finger path received signal.
        H2_est = H(:,:,2);
        G2 = (H2_est'*H2_est + (1/linearSnr)*eye(K))^-1;
        G2 = G2*H2_est';         
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mmse_le = (r_eq1 + r_eq2)/2;       
        
        % Iterate over all Tx antennas.
        for jj=1:1:K
            demoded_mmse_le = demodulate(hDemod, r_eq_sum_mmse_le(jj,:));
            E_mmse_le(:,jj:K:end) = reshape(demoded_mmse_le, modOrd, numSym/K);
        end
        % *****************************************************************************
        
        %% ********* 4) Matched Filter Bound (MFB) or Single User detection **********
        
        % Apply MRC to first finger path received signal.
        G1 = H(:,1,1)';
        r_eq1_mfb = G1*ofdm_rx_finger1_mfb;
        
        % Apply MRC to second finger path received signal.
        G2 = H(:,1,2)';
        r_eq2_mfb = G2*ofdm_rx_finger2_mfb;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mfb = (r_eq1_mfb + r_eq2_mfb)/2;        
        
        % Demodulate signal.
        E_mfb = demodulate(hDemod, r_eq_sum_mfb);
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
fileName = sprintf('M_MU_MIMO_M%s_K%s_multFad_woLargeScl_rake_perfect_est_%s.fig',m_rx_antennas,m_tx_antennas,timeStamp);
savefig(h,fileName);

% Save workspace to MAT-file.
clear h
fileName = sprintf('M_MU_MIMO_M%s_K%s_multFad_woLargeScl_rake_perfect_est_%s.mat',m_rx_antennas,m_tx_antennas,timeStamp);
save(fileName);
