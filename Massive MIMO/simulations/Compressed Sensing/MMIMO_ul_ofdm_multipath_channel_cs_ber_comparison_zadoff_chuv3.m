clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 37;                                                                                           % Number of pilots per OFDM symbol. (use prime numbers)

L = 1;                                                                                              % Number of cells.
K = 10;                                                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*(NFFT - K*Np);                                                                           % Number of symbols, i.e., number of terminals. The number of pilots must be leaved out.
NCP = 512;                                                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = true;                                                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 0:2:10;                                                                                   % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10((NFFT-K*Np+Np)/NFFT) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                                                  % Radius given in meters.
cellHole = 100;                                                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e7;
nErrors = 1000000;
debug = false;

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 0.977]*1e-6;                                             % Delay in microseconds.
gain  = [-3.010299956639812 -3.010299956639812];                    % Gain in dB (indoor).
numPaths = length(delay);                                           % Number of paths per channel.
totalNumPaths = M*K*numPaths;                                       % Total number of paths between the various antennas and Base Stations.

%% ----------------------------------- Setup Channel -----------------------------------
% Vetor de ganhos.
pos = round(delay/Ts)+1;                                            % Effective position of taps within the vector.
g = zeros(1, round(delay(end)/Ts)+1);                               % +1 is used to include the delay 0 into the vector.
for n = 1:length(delay)
    g( pos(n) ) = sqrt(10^( gain(n)/10 ));
end
delaySpreadMax = pos(length(pos));

fc = 1e9;                                                           % Carrier Freq. in MHz.
c = 3e8;                                                            % Light speed in m/s.
v = 30;                                                             % Speed in m/s.
Fs = 30.72e6;                                                       % Sampling Freq. in MHz.
Ts = 1/Fs;                                                          % Sampling period in seconds.

fd = (v*fc)/c;                                                      % Doppler frequency in Hz.

Pd = 0;                                                             % Relative power in dB.

% Parameters used for generating the Mulipath Masive MIMO Channel.
Fs_chann = 500;                                                     % Channel Sampling Rate in Hz. The channel is sampled every 2 ms. (Periodo de amostragem do canal ~1 amostra / 2*subframes)
Ts_chann = 1/Fs_chann;
N_chann = 256;                                                      % Number of samples used to sample the channel. Duration of the channel in  number of samples.
delta_f_chann = Fs_chann/N_chann;                                   % in Hz.
f = -Fs_chann/2:delta_f_chann:Fs_chann/2;
f_idx = find(f<fd);
f = f(f_idx);
f_idx = find(f>-fd);
f = f(f_idx);
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

%% @@@@@@@@@@@@@@@@ Create Pilots @@@@@@@@@@@@@@@@@@@@@@@@@@@@
u = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571];
n = (0:1:Np-1);
P = complex(zeros(K,Np),zeros(K,Np));
for u_idx=1:1:K
    P(u_idx,:) = exp((-1i.*pi.*u(u_idx).*n.*(n+1))./Np);              % Normalized pilot signal, i.e., unit power.
end

F = fft(eye(NFFT));
F = F(:,1:NCP);

%------------ Retrieve Pilot Positions ------------
if(Np==16) % Funciona somente com OMP.
    [ppos, flag_data] = getOptimumPpos16();
elseif(Np==17) % Funciona somente com OMP.
    [ppos, flag_data] = getOptimumPpos17();
elseif(Np==31) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos31();
elseif(Np==32) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos32();
elseif(Np==34) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos34();
elseif(Np==35) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos35();
elseif(Np==37) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos37();
elseif(Np==36) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos36();
elseif(Np==40) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos40();
elseif(Np==50) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos50();
elseif(Np==53) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos53();
elseif(Np==60) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos60();
elseif(Np==73) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos73();
elseif(Np==80) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos80();
elseif(Np==100) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos100();
elseif(Np==101) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos101();
elseif(Np==151) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos151();
else
    error('Invalid value for Np!!!!!');
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
[BER_MMSE_LS, BER_MMSE_MMSE, BER_MMSE_OMP, BER_MMSE_OPT] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
title(strTitle);

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
est_error_omp  = zeros(1,length(snr));
est_error_ls   = zeros(1,length(snr));
est_error_mmse = zeros(1,length(snr));
avg_error_omp  = zeros(1,length(snr));
avg_error_ls   = zeros(1,length(snr));
avg_error_mmse = zeros(1,length(snr));
N_chann = 10;
ChannelOrder = pos(length(pos))-1;
for snr_idx = 1:1:length(snr)
    
    linearSnr = 10^(snr(snr_idx)/10);
    
    ofdm_symbol_number = 0;
    idx_ch = 1;
    iter = 0;
    nBits = 0;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    nErrs_mmse_ls = 0;
    nErrs_mmse_mmse = 0;
    nErrs_mmse_omp = 0;
    nErrs_mmse_opt = 0;
    
    while(((nErrs_mmse_ls < nErrors) || (nErrs_mmse_mmse < nErrors) || (nErrs_mmse_omp < nErrors) || (nErrs_mmse_opt < nErrors)) && (nBits < nTotalOfBits))
        
        ofdm_symbol_number = ofdm_symbol_number + 1;
        
        iter = iter + 1;
        
        %---------- Transmission (UE) ----------
        % Create array of bits to modulate.
        msg = randi(hStr, [0 1], modOrd, numSym);
        
        % Modulate data.
        source = modulate(hMod, msg);
        
        % Split source among K terminals.
        Tx = complex(zeros(K,NFFT),zeros(K,NFFT));
        Tx(~flag_data) = reshape(source, K, numel(source)/K); clear source;
        for l_idx=1:1:K
            Tx(l_idx,ppos(l_idx,:)) = P(l_idx,:);
        end
        
        % Create OFDM symbol.
        sequence = (NFFT/sqrt(NFFT-(Np*K)+Np))*ifft(Tx,NFFT,2);
        
        % Add CP.
        ofdm = [sequence(:,NFFT-NCP+1:end), sequence];
        
        %---------- Multipath Channel plus Noise ----------
        H = reshape(ch(idx_ch,:), M, K, length(pos));
        
        x = [H(:,:,1)*ofdm complex(zeros(M,(pos(length(pos))-1)),zeros(M,(pos(length(pos))-1)))];
        for k = 2:length(pos)
            aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*ofdm];
            x = x + aux;
        end
        
        if(add_awgn)
            % Add channel noise power to faded data.
            r = awgn(x, snr(snr_idx), 0, hStr);
        else
            r = x;
        end
        
        %------------------- Reception (base station) ---------------------
        % Get a whole subframe.
        rx = r(:,1:end-(pos(length(pos))-1));
        
        % Remove CP.
        rx = rx(:,NCP+1:end);
        
        % Retrieve modulation symbols by applying FFT to received signal.
        ofdm_rx = (sqrt(NFFT-(Np*K)+Np)/NFFT)*fft(rx,NFFT,2);
        
        % ------------- Massive MIMO Channel Estimation -------------------
        if(ofdm_symbol_number==1)
            H_hat_omp = zeros(M, K, length(pos));
            H_hat_ls = zeros(M, K, length(pos));
            H_hat_mmse = zeros(M, K, length(pos));
            for m_idx=1:1:M
                for k_idx=1:1:K
                    
                    Pdiag = diag(P(k_idx,:));
                    
                    %% ********* Fourier Basis *********
                    Fl = F(ppos(k_idx,:).',:);
                    
                    y = ofdm_rx(m_idx,ppos(k_idx,:)).';
                    
                    % --------- Iterative Greedy Algorithm ---------
                    A = Pdiag*Fl;
                    
                    % OMP.
                    g_hat_omp = OMP_origv1(A,y,numPaths);
                    
                    % --------- Linear Algorithms ---------
                    A = Pdiag*Fl(:,1:pos(length(pos)));
                    
                    % Least Squares (LS).
                    invMatLs = (((A'*A)^(-1))*A');
                    g_hat_ls = invMatLs*y;
                    
                    % MMSE.
                    invMatMmse = (A'*A + (numPaths/linearSnr)*eye(delaySpreadMax))^-1;
                    invMatMmse = invMatMmse*A';
                    g_hat_mmse = invMatMmse*y;
                    
                    % --------- Estimate of H ---------
                    H_hat_omp(m_idx,k_idx,:)  = g_hat_omp(pos);
                    H_hat_ls(m_idx,k_idx,:)   = g_hat_ls(pos);
                    H_hat_mmse(m_idx,k_idx,:) = g_hat_mmse(pos);
                end
            end
        end
        %------------------------------------------------------------------
        
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
        ofdm_rx_finger1 = (sqrt(NFFT-(Np*K)+Np)/NFFT)*fft(rx_finger1,NFFT,2);
        
        % Retrieve modulation symbols from second finger path.
        ofdm_rx_finger2 = (sqrt(NFFT-(Np*K)+Np)/NFFT)*fft(rx_finger2,NFFT,2);
        
        %% ****************************** 1) MMSE OPT receiver *************************
        
        % Apply MMSE-LE to fisrt finger path received signal.
        H1_est = H(:,:,1);
        G1 = (((H1_est'*H1_est) + (1/linearSnr)*eye(K))^(-1));
        G1 = G1*H1_est';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MMSE-LE to second finger path received signal.
        H2_est = H(:,:,2);
        G2 = (((H2_est'*H2_est) + (1/linearSnr)*eye(K))^(-1));
        G2 = G2*H2_est';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mmse_opt = (r_eq1 + r_eq2)/2;       
        
        % Retrieve modulation symbols.
        data_symbols_mmse_opt = r_eq_sum_mmse_opt(~flag_data).';
        
        % Demodulate signal.
        E_mmse_opt = demodulate(hDemod, data_symbols_mmse_opt);
        % *****************************************************************************        
        
        %% ****************************** 2) MMSE LS receiver *************************
        
        % Apply MMSE-LE to fisrt finger path received signal.
        H1_est = H_hat_ls(:,:,1);
        G1 = (((H1_est'*H1_est) + (1/linearSnr)*eye(K))^(-1));
        G1 = G1*H1_est';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MMSE-LE to second finger path received signal.
        H2_est = H_hat_ls(:,:,2);
        G2 = (((H2_est'*H2_est) + (1/linearSnr)*eye(K))^(-1));
        G2 = G2*H2_est';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mmse_ls = (r_eq1 + r_eq2)/2;
        
        % Retrieve modulation symbols.
        data_symbols_mmse_ls = r_eq_sum_mmse_ls(~flag_data).';
        
        % Demodulate signal.
        E_mmse_ls = demodulate(hDemod, data_symbols_mmse_ls);
        % *****************************************************************************
        
        %% ****************************** 3) MMSE MMSE receiver *************************
        
        % Apply MMSE-LE to fisrt finger path received signal.
        H1_est = H_hat_mmse(:,:,1);
        G1 = (((H1_est'*H1_est) + (1/linearSnr)*eye(K))^(-1));
        G1 = G1*H1_est';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MMSE-LE to second finger path received signal.
        H2_est = H_hat_mmse(:,:,2);
        G2 = (((H2_est'*H2_est) + (1/linearSnr)*eye(K))^(-1));
        G2 = G2*H2_est';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mmse_mmse = (r_eq1 + r_eq2)/2;
              
        % Retrieve modulation symbols.
        data_symbols_mmse_mmse = r_eq_sum_mmse_mmse(~flag_data).';
        
        % Demodulate signal.
        E_mmse_mmse = demodulate(hDemod, data_symbols_mmse_mmse);
        % *****************************************************************************
        
        %% ****************************** 4) MMSE OMP receiver *************************
        
        % Apply MMSE-LE to fisrt finger path received signal.
        H1_est = H_hat_omp(:,:,1);
        G1 = (((H1_est'*H1_est) + (1/linearSnr)*eye(K))^(-1));
        G1 = G1*H1_est';
        r_eq1 = G1*ofdm_rx_finger1;
        
        % Apply MMSE-LE to second finger path received signal.
        H2_est = H_hat_omp(:,:,2);
        G2 = (((H2_est'*H2_est) + (1/linearSnr)*eye(K))^(-1));
        G2 = G2*H2_est';
        r_eq2 = G2*ofdm_rx_finger2;
        
        % Combine all the path signals by averaging them.
        r_eq_sum_mmse_omp = (r_eq1 + r_eq2)/2;
               
        % Retrieve modulation symbols.
        data_symbols_mmse_omp = r_eq_sum_mmse_omp(~flag_data).';
        
        % Demodulate signal.
        E_mmse_omp = demodulate(hDemod, data_symbols_mmse_omp);
        % *****************************************************************************
        
        %% -------------------------------- Collect errors ----------------------------
        nErrs_mmse_ls   = nErrs_mmse_ls + biterr(msg, E_mmse_ls);
        nErrs_mmse_mmse = nErrs_mmse_mmse + biterr(msg, E_mmse_mmse);
        nErrs_mmse_omp  = nErrs_mmse_omp + biterr(msg, E_mmse_omp);
        nErrs_mmse_opt  = nErrs_mmse_opt + biterr(msg, E_mmse_opt);
        
        nBits = nBits + length(msg(:));
        fprintf(1,'SNR: %d\n',snr(snr_idx));
        fprintf(1,'BER MMSE LS:   %f - nErrs_mmse_ls:   %d - nBits: %d - iter: %d\n',(nErrs_mmse_ls./nBits),nErrs_mmse_ls,nBits,iter);
        fprintf(1,'BER MMSE MMSE: %f - nErrs_mmse_mmse: %d - nBits: %d - iter: %d\n',(nErrs_mmse_mmse./nBits),nErrs_mmse_mmse,nBits,iter);
        fprintf(1,'BER MMSE OMP:  %f - nErrs_mmse_omp:  %d - nBits: %d - iter: %d\n',(nErrs_mmse_omp./nBits),nErrs_mmse_omp,nBits,iter);
        fprintf(1,'BER MMSE OPT:  %f - nErrs_mmse_opt:  %d - nBits: %d - iter: %d\n',(nErrs_mmse_opt./nBits),nErrs_mmse_opt,nBits,iter);
        fprintf(1,'\n');
        
        %% --------- Change multipath channel according to its sampling rate. ---------
        if(ofdm_symbol_number==2*numSymbInSubframe)
            ofdm_symbol_number = 0;
            idx_ch = idx_ch + 1;
            if(idx_ch > N_chann)
                idx_ch = 1;
            end
        end
        %------------------------------------------------------------------
        
    end
    
    % Calculate BER for current SNR point.
    BER_MMSE_LS(snr_idx) = nErrs_mmse_ls./nBits;
    BER_MMSE_MMSE(snr_idx) = nErrs_mmse_mmse./nBits;
    BER_MMSE_OMP(snr_idx) = nErrs_mmse_omp./nBits;
    BER_MMSE_OPT(snr_idx) = nErrs_mmse_opt./nBits;
    
    % Plot results.
    semilogy(EbNoVec(1:snr_idx), BER_MMSE_LS(1:snr_idx), 'ks', ...
        EbNoVec(1:snr_idx), BER_MMSE_MMSE(1:snr_idx), 'b*', ...
        EbNoVec(1:snr_idx), BER_MMSE_OMP(1:snr_idx), 'ro', ...
        EbNoVec(1:snr_idx), BER_MMSE_OPT(1:snr_idx), 'k');
    legend('MMSE LS', 'MMSE MMSE', 'MMSE OMP', 'MMSE OPT');
    drawnow;
    
end

% Draw the lines.
semilogy(EbNoVec, BER_MMSE_LS, 'k-', EbNoVec, BER_MMSE_MMSE, 'b-', EbNoVec, BER_MMSE_OMP, 'r-', EbNoVec, BER_MMSE_OPT, 'k-');
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
    
    n_pilots = '';
    for j=1:length(Np)
        n_pilots = strcat(n_pilots, sprintf('_%d',Np(j)));
    end
    
    % Get timestamp for saving files.
    timeStamp = datestr(now,30);
    
    % Save figure to FIG-file.
    fileName = sprintf('Massive_MU_MIMO_M%s_K%s_Np%s_BER_comparison_%s.fig',m_rx_antennas,m_tx_antennas,n_pilots,timeStamp);
    savefig(figura,fileName); clear figura
    
    % Save workspace to MAT-file.
    fileName = sprintf('Massive_MU_MIMO_M%s_K%s_Np%s_BER_comparison_%s.mat',m_rx_antennas,m_tx_antennas,n_pilots,timeStamp);
    save(fileName);
end