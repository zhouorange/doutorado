clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 37;                                                                                            % Number of pilots per OFDM symbol. (use prime numbers)

L = 1;                                                                                              % Number of cells.
K = 10;                                                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*NFFT;                                                                                    % Number of symbols, i.e., number of terminals. The number of pilots must be left out.
NCP = 512;                                                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = false;                                                                                   % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 10000; %-20; %-10:2:6;%-20:2:40;                                                          % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);                                    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                                                  % Radius given in meters.
cellHole = 100;                                                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e9;
nErrors = 1e8;
debug = false;

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 30]*Ts;                                                  % Delay in microseconds.
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
P = complex(zeros(K,Np),zeros(K,Np));
for u_idx=1:1:K
    % BPSK pilots.
    ip = rand(1,Np) > 0.5;
    P(u_idx,:) = 2*ip-1;
end

F = fft(eye(NFFT));
F = F(:,1:NCP);

%------------ Retrieve Best Pilot Positions ------------
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

% Pre-allocate variables for speed.
dist = zeros(length(bits), 1);
[BER_LS, BER_OMP, BER_FULL, BER_ELS, BER_GAELS] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-10 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
title(strTitle);

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for snr_idx = 1:1:length(snr)
    
    linearSnr = 10^(snr(snr_idx)/10);
    
    ofdm_symbol_number = 0;
    idx_ch = 1;
    iter = 0;
    nBits = 0;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    nErrs_omp = 0;
    nErrs_ls = 0;
    nErrs_gaels = 0;
    nErrs_els = 0;
    nErrs_full = 0;
    
    while(((nErrs_ls < nErrors) || (nErrs_omp < nErrors) || (nErrs_full < nErrors) || (nErrs_els < nErrors) || (nErrs_gaels < nErrors)) && (nBits < nTotalOfBits))
        
        ofdm_symbol_number = ofdm_symbol_number + 1;
        
        iter = iter + 1;
        
        %---------- Transmission (UE) ----------
        Tx = complex(zeros(K,NFFT),zeros(K,NFFT));
        
        if(ofdm_symbol_number==1)
            % TX Pilots.
            for l_idx=1:1:K
                Tx(l_idx,ppos(l_idx,:)) = P(l_idx,:);
            end
        else
            % TX Data.
            % Create array of bits to modulate.
            msg = randi(hStr, [0 1], modOrd, numSym);
            
            % Modulate data.
            source = modulate(hMod, msg);
            
            % Split source among K terminals.
            Tx = reshape(source, K, numel(source)/K); clear source;
        end
        
        % Create OFDM symbol.
        sequence = (NFFT/sqrt(NFFT))*ifft(Tx,NFFT,2);
        
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
        ofdm_rx = (sqrt(NFFT)/NFFT)*fft(rx,NFFT,2);
        
        % **************** Massive MIMO Channel Estimation ****************
        if(ofdm_symbol_number==1)
            error_omp = 0;
            error_ls = 0;
            error_gaels = 0;
            error_els = 0;
            H_hat_omp = zeros(M, K, length(pos));
            H_hat_ls = zeros(M, K, length(pos));
            H_hat_gaels = zeros(M, K, length(pos));
            H_hat_els = zeros(M, K, length(pos));
            H_hat_fd_omp = zeros(M, K, NFFT);
            H_hat_fd_ls = zeros(M, K, NFFT);
            H_hat_fd_gaels = zeros(M, K, NFFT);
            H_hat_fd_els = zeros(M, K, NFFT);
            H_full_fd = zeros(M, K, NFFT);
            for m_idx=1:1:M
                for k_idx=1:1:K
                    
                    Pdiag = diag(P(k_idx,:));
                    
                    %% ********* Fourier Basis *********
                    Fl = F(ppos(k_idx,:).',:);
                    
                    y = ofdm_rx(m_idx,ppos(k_idx,:)).';
                    
                    % -- Iterative Compressive Sensing Greedy Algorithm ---
                    A = Pdiag*Fl;
                    
                    % OMP.
                    g_hat_omp = OMP_origv1(A,y,numPaths);
                    
                    % --------- Linear Algorithms ---------
                    A = Pdiag*Fl(:,1:pos(length(pos)));
                    
                    % Least Squares (LS).
                    invMatLs = (((A'*A)^(-1))*A');
                    g_hat_ls = invMatLs*y;
                    
                    % ELS: Enhanced Least Squares.
                    % Forma que melhora a performance do LS, porém, deve-se saber a posição dos taps diferentes de zero.
                    % A posição é encontrada através do Matching Pursuit.
                    [g_hat_mp tap_pos tap_res]= mp(A,y); % Matching Pursuit.
                    tap_pos = sort(tap_pos(1:2));
                    d = zeros(1,pos(length(pos)));
                    d(tap_pos) = [1 1];
                    D = diag(d);
                    B = Pdiag*Fl(:,1:pos(length(pos)))*D;
                    Bl = B(:,[tap_pos(1) tap_pos(length(tap_pos))]);
                    invMatEls = (((Bl'*Bl)^(-1))*Bl');
                    g_hat_els_aux = invMatEls*y;
                    g_hat_els = complex(zeros(pos(length(pos)),1),zeros(pos(length(pos)),1));
                    g_hat_els(tap_pos) = g_hat_els_aux;
                    
                    % GA-ELS: Genie Aided Enhanced Least Squares.
                    % Forma que melhora a performance do MMSE, porém, deve-se saber a posição dos taps diferentes de zero.
                    d = zeros(1,pos(length(pos)));
                    d(pos) = [1 1];
                    D = diag(d);
                    B = Pdiag*Fl(:,1:pos(length(pos)))*D;
                    Bl = B(:,[pos(1) pos(length(pos))]);
                    invMatGAEls = (((Bl'*Bl)^(-1))*Bl');
                    g_hat_gaels = invMatGAEls*y;
                    g_hat_gaels = [g_hat_gaels(1); zeros(31-2,1); g_hat_gaels(2)];
                    
                    % --------- Estimate of H (Frequency Domain) ---------
                    H1Ant_omp = fft([g_hat_omp; zeros(NFFT-length(g_hat_omp),1)],NFFT);
                    H1Ant_ls = fft([g_hat_ls; zeros(NFFT-length(g_hat_ls),1)],NFFT);
                    H1Ant_gaels = fft([g_hat_gaels; zeros(NFFT-length(g_hat_gaels),1)],NFFT);
                    H1Ant_els = fft([g_hat_els; zeros(NFFT-length(g_hat_els),1)],NFFT);
                    
                    h_orig = zeros(NFFT,1);
                    h_orig(pos,1) = [H(m_idx,k_idx,1); H(m_idx,k_idx,2)];
                    H_orig = fft(h_orig,NFFT);
                    
                    H_hat_fd_omp(m_idx,k_idx,:)   = H1Ant_omp;
                    H_hat_fd_ls(m_idx,k_idx,:)    = H1Ant_ls;
                    H_hat_fd_gaels(m_idx,k_idx,:) = H1Ant_gaels;
                    H_hat_fd_els(m_idx,k_idx,:)   = H1Ant_els;
                    H_full_fd(m_idx,k_idx,:)      = H_orig;
                    
                    % --------- Estimate of H (Time Domain) ---------
                    H_hat_omp(m_idx,k_idx,:)   = g_hat_omp(pos);
                    H_hat_ls(m_idx,k_idx,:)    = g_hat_ls(pos);
                    H_hat_gaels(m_idx,k_idx,:) = g_hat_gaels(pos);
                    H_hat_els(m_idx,k_idx,:)   = g_hat_els(pos);
                    
                    % ---------- Calculate Freq. Domain Channel Estimation Error -----------
                    error_omp = error_omp + (sum(abs(H_orig-H1Ant_omp).^2)/NFFT);
                    error_ls = error_ls + (sum(abs(H_orig-H1Ant_ls).^2)/NFFT);
                    error_gaels = error_gaels + (sum(abs(H_orig-H1Ant_gaels).^2)/NFFT);
                    error_els = error_els + (sum(abs(H_orig-H1Ant_els).^2)/NFFT);
                end
            end
            % --------- Print Frequency Domain Estimation Error -----------
            error_omp = error_omp/(M*K);
            error_ls = error_ls/(M*K);
            error_gaels = error_gaels/(M*K);
            error_els = error_els/(M*K);
            fprintf(1,'OMP FD Error Estimation    : %d\n',error_omp);
            fprintf(1,'LS FD Error Estimation     : %d\n',error_ls);
            fprintf(1,'GA-ELS FD Error Estimation : %d\n',error_gaels);
            fprintf(1,'MP-ELS FD Error Estimation : %d\n',error_els);
            fprintf(1,'\n');            
        else           
            % ****** Apply Frequency Domain Equalization per subcarrier *******
            rx_omp = zeros(K,NFFT);
            rx_ls = zeros(K,NFFT);
            rx_gaels = zeros(K,NFFT);
            rx_els = zeros(K,NFFT);
            rx_full = zeros(K,NFFT);
            for sc=1:1:NFFT
                % ----------- 1) Estimation: OMP - Equalizer: LS --------------
                H_omp = H_hat_fd_omp(:,:,sc);
                H_omp_inv = (((H_omp'*H_omp)^(-1))*H_omp');
                rx_omp(:,sc) = H_omp_inv*ofdm_rx(:,sc);
                
                % ----------- 2) Estimation: LS - Equalizer: LS ---------------
                H_ls = H_hat_fd_ls(:,:,sc);
                H_ls_inv = (((H_ls'*H_ls)^(-1))*H_ls');
                rx_ls(:,sc) = H_ls_inv*ofdm_rx(:,sc);
                
                % ----------- 3) Estimation: GA-MP-LS - Equalizer: LS ------------
                H_gaels = H_hat_fd_gaels(:,:,sc);
                H_gaels_inv = (((H_gaels'*H_gaels)^(-1))*H_gaels');
                rx_gaels(:,sc) = H_gaels_inv*ofdm_rx(:,sc);
                
                % ----------- 4) Estimation: MP-LS - Equalizer: LS ------------
                H_els = H_hat_fd_els(:,:,sc);
                H_els_inv = (((H_els'*H_els)^(-1))*H_els');
                rx_els(:,sc) = H_els_inv*ofdm_rx(:,sc);
                
                % ----------- 5) Estimation: FULL - Equalizer: LS -------------
                H_fd = H_full_fd(:,:,sc);
                H_fd_inv = (((H_fd'*H_fd)^(-1))*H_fd');
                rx_full(:,sc) = H_fd_inv*ofdm_rx(:,sc);
            end
            
            % *********************** Demodulate signal ***********************
            % -------------- 1) Estimation: OMP - Equalizer: LS ---------------
            data_symbols_omp = reshape(rx_omp,1,NFFT*K);
            E_omp = demodulate(hDemod, data_symbols_omp);
            
            % -------------- 2) Estimation: LS - Equalizer: LS ----------------
            data_symbols_ls = reshape(rx_ls,1,NFFT*K);
            E_ls = demodulate(hDemod, data_symbols_ls);
            
            % -------------- 3) Estimation: GA-LS - Equalizer: LS -------------
            data_symbols_gaels = reshape(rx_gaels,1,NFFT*K);
            E_gaels = demodulate(hDemod, data_symbols_gaels);
            
            % -------------- 4) Estimation: MP-LS - Equalizer: LS -------------
            data_symbols_els = reshape(rx_els,1,NFFT*K);
            E_els = demodulate(hDemod, data_symbols_els);
            
            % -------------- 5) Estimation: Full - Equalizer: LS --------------
            data_symbols_full = reshape(rx_full,1,NFFT*K);
            E_full = demodulate(hDemod, data_symbols_full);
            
            %% -------------------------------- Collect errors ----------------------------
            nErrs_omp = nErrs_omp + biterr(msg, E_omp);
            nErrs_ls = nErrs_ls + biterr(msg, E_ls);
            nErrs_gaels = nErrs_gaels + biterr(msg, E_gaels);
            nErrs_els = nErrs_els + biterr(msg, E_els);
            nErrs_full = nErrs_full + biterr(msg, E_full);
            
            nBits = nBits + length(msg(:));
            %         fprintf(1,'SNR: %d\n',snr(snr_idx));
            %         fprintf(1,'BER LS-LS:     %f - nErrs_ls:    %d - nBits: %d - iter: %d\n',(nErrs_ls./nBits),nErrs_ls,nBits,iter);
            %         fprintf(1,'BER GA-ELS-LS: %f - nErrs_gaels: %d - nBits: %d - iter: %d\n',(nErrs_gaels./nBits),nErrs_gaels,nBits,iter);
            %         fprintf(1,'BER MP-LS-LS:  %f - nErrs_els:   %d - nBits: %d - iter: %d\n',(nErrs_els./nBits),nErrs_els,nBits,iter);
            %         fprintf(1,'BER OMP-LS:    %f - nErrs_omp:   %d - nBits: %d - iter: %d\n',(nErrs_omp./nBits),nErrs_omp,nBits,iter);
            %         fprintf(1,'BER FULL-LS:   %f - nErrs_full:  %d - nBits: %d - iter: %d\n',(nErrs_full./nBits),nErrs_full,nBits,iter);
            %         fprintf(1,'\n');
            
        end
        
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
    BER_LS(snr_idx) = nErrs_ls./nBits;
    BER_OMP(snr_idx) = nErrs_omp./nBits;
    BER_FULL(snr_idx) = nErrs_full./nBits;
    BER_ELS(snr_idx) = nErrs_els./nBits;
    BER_GAELS(snr_idx) = nErrs_gaels./nBits;
    
    % Plot results.
    semilogy(EbNoVec(1:snr_idx), BER_LS(1:snr_idx), 'ks', ...
        EbNoVec(1:snr_idx), BER_GAELS(1:snr_idx), 'k*', ...
        EbNoVec(1:snr_idx), BER_ELS(1:snr_idx), 'ko', ...
        EbNoVec(1:snr_idx), BER_OMP(1:snr_idx), 'ro', ...
        EbNoVec(1:snr_idx), BER_FULL(1:snr_idx), 'g');
    legend('LS-LS', 'GA-ELS-LS', 'MP-LS', 'OMP-LS', 'FULL-LS');
    drawnow;
    
end

% Draw the lines.
semilogy(EbNoVec, BER_LS, 'k-', EbNoVec, BER_GAELS, 'k-', EbNoVec, BER_ELS, 'k-', EbNoVec, BER_OMP, 'r-', EbNoVec, BER_FULL, 'g-');
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
    fileName = sprintf('Massive_MU_MIMO_M%s_K%s_Np%s_BER_comparison_%s_mp_ls_v6.fig',m_rx_antennas,m_tx_antennas,n_pilots,timeStamp);
    savefig(figura,fileName); clear figura
    
    % Save workspace to MAT-file.
    fileName = sprintf('Massive_MU_MIMO_M%s_K%s_Np%s_BER_comparison_%s_mp_ls_v6.mat',m_rx_antennas,m_tx_antennas,n_pilots,timeStamp);
    save(fileName);
end