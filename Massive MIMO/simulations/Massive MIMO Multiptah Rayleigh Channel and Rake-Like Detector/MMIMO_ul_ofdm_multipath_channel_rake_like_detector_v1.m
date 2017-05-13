clear all;close all;clc

script_version = 1;

MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID);

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

pilot_type = 'bpsk';                                                                                % Defines the type of pilot to be used: bpsk or zadoffchu

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
numSymbInSubframe = 50;                                                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = true;                                                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = -20:2:20;                                                                                  % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);                                    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                                                  % Radius given in meters.
cellHole = 100;                                                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                                                          % enable/disable power control at eNodeB.

% %               -20 -18 -16 -14 -12 -10 -8  -6  -4  -2  0
% nTotalOfBits = [1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e8 1e8 1e8 1e8];
% nErrors =      [1e6 1e6 1e6 1e6 1e6 1e6 1e6 1e7 1e7 1e7 1e7];

%               -20 -18 -16 -14 -12 -10 -8  -6  -4  -2  0   2   4   6   8   10  12  14  16  18  20
nTotalOfBits = [1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8 1e8];
nErrors =      [1e6 1e6 1e6 1e6 1e6 1e6 1e6 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7 1e7];

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
R_hh = diag(g.^2);
h_var = g(1)^2;
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
if(strcmp(pilot_type,'bpsk'))
    P = complex(zeros(K,Np),zeros(K,Np));
    for u_idx=1:1:K
        % BPSK pilots.
        ip = rand(1,Np) > 0.5;
        P(u_idx,:) = 2*ip-1;
    end
elseif(strcmp(pilot_type,'zadoffchu'))
    u = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571];
    n = (0:1:Np-1);
    P = complex(zeros(K,Np),zeros(K,Np));
    for u_idx=1:1:K
        P(u_idx,:) = exp((-1i.*pi.*u(u_idx).*n.*(n+1))./Np);              % Normalized pilot signal, i.e., unit power.
    end
else
    error('Undefined pilot type');
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
[BER_LS, BER_MMSE, BER_OMP, BER_FULL, BER_ELS, BER_GAELS, BER_COSAMP, BER_MP] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-10 1]);
xlabel('Eb/No [dB]'); ylabel('BER'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Uplink MU-MMIMO Rake-Like Detector - M: %d - K: %d - Np: %d',M,K,Np);
title(strTitle);

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
error_td_omp = zeros(1,length(snr));
error_td_cosamp = zeros(1,length(snr));
error_td_mp = zeros(1,length(snr));
error_td_ls = zeros(1,length(snr));
error_td_gaels = zeros(1,length(snr));
error_td_els = zeros(1,length(snr));
error_td_mmse = zeros(1,length(snr));
a = (NFFT/sqrt(NFFT));
Es = 1; % BPSK symbol energy.
for snr_idx = 1:1:length(snr)
    
    linearSnr = 10^(snr(snr_idx)/10);
    
    matrixFactor = ((NFFT/linearSnr)/(a*a*Es))*eye(K);
    
    estimation_counter = 0;
    ofdm_symbol_number = 0;
    idx_ch = 1;
    iter = 0;
    nBits = 0;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    nErrs_omp_mmse = 0;
    nErrs_cosamp_mmse = 0;
    nErrs_mp_mmse = 0;
    nErrs_ls_mmse = 0;
    nErrs_mmse_mmse = 0;
    nErrs_gaels_mmse = 0;
    nErrs_els_mmse = 0;
    nErrs_full_mmse = 0;
    
    while(((nErrs_ls_mmse < nErrors(snr_idx)) || (nErrs_omp_mmse < nErrors(snr_idx)) || (nErrs_full_mmse < nErrors(snr_idx)) || (nErrs_els_mmse < nErrors(snr_idx)) || (nErrs_gaels_mmse < nErrors(snr_idx)) || (nErrs_mmse_mmse < nErrors(snr_idx)) || (nErrs_cosamp_mmse < nErrors(snr_idx)) || (nErrs_mp_mmse < nErrors(snr_idx))) && (nBits < nTotalOfBits(snr_idx)))
        
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
            error_cosamp = 0;
            error_mp = 0;
            error_ls = 0;
            error_gaels = 0;
            error_els = 0;
            error_mmse = 0;
            H_hat_omp = zeros(M, K, length(pos));
            H_hat_cosamp = zeros(M, K, length(pos));
            H_hat_mp = zeros(M, K, length(pos));
            H_hat_ls = zeros(M, K, length(pos));
            H_hat_gaels = zeros(M, K, length(pos));
            H_hat_els = zeros(M, K, length(pos));
            H_hat_mmse = zeros(M, K, length(pos));
            H_hat_fd_omp = zeros(M, K, NFFT);
            H_hat_fd_cosamp = zeros(M, K, NFFT);
            H_hat_fd_mp = zeros(M, K, NFFT);
            H_hat_fd_ls = zeros(M, K, NFFT);
            H_hat_fd_gaels = zeros(M, K, NFFT);
            H_hat_fd_els = zeros(M, K, NFFT);
            H_hat_fd_mmse = zeros(M, K, NFFT);
            H_full_fd = zeros(M, K, NFFT);
            for m_idx=1:1:M
                for k_idx=1:1:K
                    
                    Pdiag = diag(P(k_idx,:));
                    
                    %% ********* Fourier Basis *********
                    Fl = F(ppos(k_idx,:).',:);
                    
                    y = ofdm_rx(m_idx,ppos(k_idx,:)).';
                    
                    % -- Iterative Compressive Sensing: Greedy Algorithms ---
                    A = Pdiag*Fl;
                    
                    % OMP.
                    g_hat_omp = OMP_origv1(A,y,numPaths); % OMP original
                    
                    % CoSamp.
                    g_hat_cosamp = CoSaMPv1(A,y,numPaths);
                    
                    % Matching Pursuit (MP).
                    g_hat_mp = mpv2(A,y,numPaths);
                    
                    % --------- Linear Algorithms ---------
                    A = Pdiag*Fl(:,1:pos(length(pos)));
                    
                    % Least Squares (LS).
                    invMatLs = (((A'*A)^(-1))*A');
                    g_hat_ls = invMatLs*y;
                    
                    % ELS: Enhanced Least Squares.
                    % Forma que melhora a performance do LS, porém, deve-se saber a posição dos taps diferentes de zero.
                    % A posição é encontrada através do Matching Pursuit.
                    [g_hat_mp2, tap_pos, tap_res]= mp(A,y); % Matching Pursuit.
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
                    
                    % Minimum Mean-Square Error (MMSE) Estimation.
                    % ******* Implementation that does not take into account the position of the taps. *******
                    % Nao existe diferenca entre as 2 formas seguintes de
                    % se executar o MMSE. Porém, a forma abaixo gera muitos
                    % erros quando se tem SNR alta.
                    %invMatMmse = ( (A')*((A*A' + ((NFFT/linearSnr)/(a*a*h_var))*eye(Np))^(-1)));
                    invMatMmse = (((A'*A + ((NFFT/linearSnr)/(a*a*h_var))*eye(pos(length(pos))))^(-1))*A');
                    g_hat_mmse = invMatMmse*y;
                    
                    %                     %*****************************************************************************************
                    %                     % ********* Implementation that DOES take into account the position of the taps. *********
                    %                     part0 = R_hh*(A');
                    %                     part1 = A*R_hh*(A');
                    %                     part2 = ((NFFT/linearSnr)/(a*a))*eye(Np);
                    %                     part3 = (part1 + part2)^(-1);
                    %                     invMatMmse = part0*part3;
                    %                     g_hat_mmse = invMatMmse*y;
                    %                     %*****************************************************************************************
                    
                    % --------- Estimate of H (Time Domain) ---------
                    H_hat_omp(m_idx,k_idx,:)    = g_hat_omp(pos);
                    H_hat_cosamp(m_idx,k_idx,:) = g_hat_cosamp(pos);
                    H_hat_mp(m_idx,k_idx,:)     = g_hat_mp(pos);
                    H_hat_ls(m_idx,k_idx,:)     = g_hat_ls(pos);
                    H_hat_gaels(m_idx,k_idx,:)  = g_hat_gaels(pos);
                    H_hat_els(m_idx,k_idx,:)    = g_hat_els(pos);
                    H_hat_mmse(m_idx,k_idx,:)   = g_hat_mmse(pos);
                    
                    % ---------- Calculate Time Domain Channel Estimation Error -----------
                    g_orig = zeros(NCP,1);
                    g_orig(pos) = H(m_idx,k_idx,:);
                    
                    error_omp = error_omp + (sum(abs(g_orig-g_hat_omp).^2)/length(g_hat_omp));
                    error_cosamp = error_cosamp + (sum(abs(g_orig-g_hat_cosamp).^2)/length(g_hat_cosamp));
                    error_mp = error_mp + (sum(abs(g_orig-g_hat_mp).^2)/length(g_hat_mp));
                    error_ls = error_ls + (sum(abs(g_orig(1:pos(length(pos)))-g_hat_ls).^2)/length(g_hat_ls));
                    error_gaels = error_gaels + (sum(abs(g_orig(1:pos(length(pos)))-g_hat_gaels).^2)/length(g_hat_gaels));
                    error_els = error_els + (sum(abs(g_orig(1:pos(length(pos)))-g_hat_els).^2)/length(g_hat_els));
                    error_mmse = error_mmse + (sum(abs(g_orig(1:pos(length(pos)))-g_hat_mmse).^2)/length(g_hat_mmse));
                end
            end
            
            % --------- Print Time Domain Estimation Error -----------
            estimation_counter = estimation_counter + 1;
            
            error_td_omp(snr_idx) = error_td_omp(snr_idx) + (error_omp/(M*K));
            error_td_cosamp(snr_idx) = error_td_cosamp(snr_idx) + (error_cosamp/(M*K));
            error_td_mp(snr_idx) = error_td_mp(snr_idx) + (error_mp/(M*K));
            error_td_ls(snr_idx) = error_td_ls(snr_idx) + (error_ls/(M*K));
            error_td_gaels(snr_idx) = error_td_gaels(snr_idx) + (error_gaels/(M*K));
            error_td_els(snr_idx) = error_td_els(snr_idx) + (error_els/(M*K));
            error_td_mmse(snr_idx) = error_td_mmse(snr_idx) + (error_mmse/(M*K));
            
            fprintf(1,'OFDM symbol number: %d\n',ofdm_symbol_number);
            fprintf(1,'OMP TD Error Estimation    : %d\n',(error_omp/(M*K)));
            fprintf(1,'CoSamp TD Error Estimation : %d\n',(error_cosamp/(M*K)));
            fprintf(1,'MP TD Error Estimation     : %d\n',(error_mp/(M*K)));
            fprintf(1,'LS TD Error Estimation     : %d\n',(error_ls/(M*K)));
            fprintf(1,'MMSE TD Error Estimation   : %d\n',(error_mmse/(M*K)));
            fprintf(1,'GA-ELS TD Error Estimation : %d\n',(error_gaels/(M*K)));
            fprintf(1,'MP-ELS TD Error Estimation : %d\n',(error_els/(M*K)));
            fprintf(1,'\n');
        else

            
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
            
            %% ****************************** 1) FULL/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H(:,:,1);
            G1 = (H1_est'*H1_est + matrixFactor )^-1;
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H(:,:,2);
            G2 = (H2_est'*H2_est + matrixFactor )^-1;
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mmse_opt = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mmse_opt = reshape(r_eq_sum_mmse_opt,1,NFFT*K);
            
            % Demodulate signal.
            E_full_mmse = demodulate(hDemod, data_symbols_mmse_opt);
            % *****************************************************************************
            
            %% ****************************** 2) LS/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_ls(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_ls(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mmse_ls = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mmse_ls = reshape(r_eq_sum_mmse_ls,1,NFFT*K);
            
            % Demodulate signal.
            E_ls_mmse = demodulate(hDemod, data_symbols_mmse_ls);
            % *****************************************************************************
            
            %% ****************************** 3) MMSE/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_mmse(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_mmse(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mmse_mmse = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mmse_mmse = reshape(r_eq_sum_mmse_mmse,1,NFFT*K);
            
            % Demodulate signal.
            E_mmse_mmse = demodulate(hDemod, data_symbols_mmse_mmse);
            % *****************************************************************************

            %% ****************************** 4) OMP/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_omp(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_omp(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mmse_omp = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mmse_omp = reshape(r_eq_sum_mmse_omp,1,NFFT*K);
            
            % Demodulate signal.
            E_omp_mmse = demodulate(hDemod, data_symbols_mmse_omp);
            % *****************************************************************************
            
            %% ****************************** 5) CoSamp/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_cosamp(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_cosamp(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mmse_cosamp = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mmse_cosamp = reshape(r_eq_sum_mmse_cosamp,1,NFFT*K);
            
            % Demodulate signal.
            E_cosamp_mmse = demodulate(hDemod, data_symbols_mmse_cosamp);
            % *****************************************************************************            
            
            %% ****************************** 6) MP/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_mp(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_mp(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_mp_mmse = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_mp_mmse = reshape(r_eq_sum_mp_mmse,1,NFFT*K);
            
            % Demodulate signal.
            E_mp_mmse = demodulate(hDemod, data_symbols_mp_mmse);
            % *****************************************************************************
            
            %% ****************************** 5) GA-ELS/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_gaels(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_gaels(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_gaels_mmse = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_gaels_mmse = reshape(r_eq_sum_gaels_mmse,1,NFFT*K);
            
            % Demodulate signal.
            E_gaels_mmse = demodulate(hDemod, data_symbols_gaels_mmse);
            % *****************************************************************************

            %% ****************************** 5) MP-LS/MMSE receiver *************************
            
            % Apply MMSE-LE to fisrt finger path received signal.
            H1_est = H_hat_els(:,:,1);
            G1 = (((H1_est'*H1_est) + matrixFactor )^(-1));
            G1 = G1*H1_est';
            r_eq1 = G1*ofdm_rx_finger1;
            
            % Apply MMSE-LE to second finger path received signal.
            H2_est = H_hat_els(:,:,2);
            G2 = (((H2_est'*H2_est) + matrixFactor )^(-1));
            G2 = G2*H2_est';
            r_eq2 = G2*ofdm_rx_finger2;
            
            % Combine all the path signals by averaging them.
            r_eq_sum_els_mmse = (r_eq1 + r_eq2)/2;
            
            % Retrieve modulation symbols.
            data_symbols_els_mmse = reshape(r_eq_sum_els_mmse,1,NFFT*K);
            
            % Demodulate signal.
            E_els_mmse = demodulate(hDemod, data_symbols_els_mmse);
            % *****************************************************************************
            
            %% -------------------------------- Collect errors ----------------------------
            nErrs_full_mmse = nErrs_full_mmse + biterr(msg, E_full_mmse);
            nErrs_ls_mmse = nErrs_ls_mmse + biterr(msg, E_ls_mmse);
            nErrs_mmse_mmse = nErrs_mmse_mmse + biterr(msg, E_mmse_mmse);
            nErrs_omp_mmse = nErrs_omp_mmse + biterr(msg, E_omp_mmse);
            nErrs_cosamp_mmse = nErrs_cosamp_mmse + biterr(msg, E_cosamp_mmse);
            nErrs_mp_mmse = nErrs_mp_mmse + biterr(msg, E_mp_mmse);
            nErrs_gaels_mmse = nErrs_gaels_mmse + biterr(msg, E_gaels_mmse);
            nErrs_els_mmse = nErrs_els_mmse + biterr(msg, E_els_mmse);
            
            nBits = nBits + length(msg(:));
            fprintf(1,'OFDM symbol number: %d\n',ofdm_symbol_number);
            fprintf(1,'SNR: %d\n',snr(snr_idx));
            fprintf(1,'Eb/N0: %d\n',EbNoVec(snr_idx));
            fprintf(1,'BER LS/MMSE:     %f - nErrs_ls:     %d - nBits: %d - iter: %d\n',(nErrs_ls_mmse./nBits),nErrs_ls_mmse,nBits,iter);
            fprintf(1,'BER MMSE/MMSE:   %f - nErrs_mmse:   %d - nBits: %d - iter: %d\n',(nErrs_mmse_mmse./nBits),nErrs_mmse_mmse,nBits,iter);
            fprintf(1,'BER GA-ELS/MMSE: %f - nErrs_gaels:  %d - nBits: %d - iter: %d\n',(nErrs_gaels_mmse./nBits),nErrs_gaels_mmse,nBits,iter);
            fprintf(1,'BER MP-LS/MMSE:  %f - nErrs_els:    %d - nBits: %d - iter: %d\n',(nErrs_els_mmse./nBits),nErrs_els_mmse,nBits,iter);
            fprintf(1,'BER OMP/MMSE:    %f - nErrs_omp:    %d - nBits: %d - iter: %d\n',(nErrs_omp_mmse./nBits),nErrs_omp_mmse,nBits,iter);
            fprintf(1,'BER COSAMP/MMSE: %f - nErrs_cosamp: %d - nBits: %d - iter: %d\n',(nErrs_cosamp_mmse./nBits),nErrs_cosamp_mmse,nBits,iter);
            fprintf(1,'BER MP/MMSE:     %f - nErrs_mp:     %d - nBits: %d - iter: %d\n',(nErrs_mp_mmse./nBits),nErrs_mp_mmse,nBits,iter);
            fprintf(1,'BER FULL/MMSE:   %f - nErrs_full:   %d - nBits: %d - iter: %d\n',(nErrs_full_mmse./nBits),nErrs_full_mmse,nBits,iter);
            fprintf(1,'\n');
            
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
    
    % Calculate Frequency Domain Estimation Error for each SNR point.
    error_td_omp(snr_idx) = error_td_omp(snr_idx)/estimation_counter;
    error_td_cosamp(snr_idx) = error_td_cosamp(snr_idx)/estimation_counter;
    error_td_mp(snr_idx) = error_td_mp(snr_idx)/estimation_counter;
    error_td_ls(snr_idx) = error_td_ls(snr_idx)/estimation_counter;
    error_td_gaels(snr_idx) = error_td_gaels(snr_idx)/estimation_counter;
    error_td_els(snr_idx) = error_td_els(snr_idx)/estimation_counter;
    error_td_mmse(snr_idx) = error_td_mmse(snr_idx)/estimation_counter;
    
    % Calculate BER for current SNR point.
    BER_LS(snr_idx) = nErrs_ls_mmse./nBits;
    BER_MMSE(snr_idx) = nErrs_mmse_mmse./nBits;
    BER_OMP(snr_idx) = nErrs_omp_mmse./nBits;
    BER_COSAMP(snr_idx) = nErrs_cosamp_mmse./nBits;
    BER_MP(snr_idx) = nErrs_mp_mmse./nBits;
    BER_FULL(snr_idx) = nErrs_full_mmse./nBits;
    BER_ELS(snr_idx) = nErrs_els_mmse./nBits;
    BER_GAELS(snr_idx) = nErrs_gaels_mmse./nBits;
    
    % Plot results.
    semilogy(EbNoVec(1:snr_idx), BER_LS(1:snr_idx), 'ks', ...
        EbNoVec(1:snr_idx), BER_MMSE(1:snr_idx), 'bs', ...
        EbNoVec(1:snr_idx), BER_GAELS(1:snr_idx), 'k*', ...
        EbNoVec(1:snr_idx), BER_ELS(1:snr_idx), 'ko', ...
        EbNoVec(1:snr_idx), BER_OMP(1:snr_idx), 'ro', ...
        EbNoVec(1:snr_idx), BER_COSAMP(1:snr_idx), 'rs', ...
        EbNoVec(1:snr_idx), BER_MP(1:snr_idx), 'gs', ...
        EbNoVec(1:snr_idx), BER_FULL(1:snr_idx), 'g');
    legend('LS/MMSE', 'MMSE/MMSE', 'GA-MP-LS/MMSE', 'MP-LS/MMSE', 'OMP/MMSE', 'COSAMP/MMSE', 'MP/MMSE', 'FULL/MMSE');
    drawnow;
    
end

% Draw the lines.
semilogy(EbNoVec, BER_LS, 'k-', EbNoVec, BER_MMSE, 'b-', EbNoVec, BER_GAELS, 'k-', EbNoVec, BER_ELS, 'k-', EbNoVec, BER_OMP, 'r-', EbNoVec, BER_COSAMP, 'r-', EbNoVec, BER_MP, 'g-', EbNoVec, BER_FULL, 'g-');
hold off;

% Plot frequency domain estimation error for each one of the techniques.
fdee_figure = figure; grid on; hold on;
xlabel('Eb/No [dB]');ylabel('MSE');
title('Time-Domain Estimation Error');

semilogy(EbNoVec, error_td_ls, 'ks', ...
    EbNoVec, error_td_mmse, 'bs', ...
    EbNoVec, error_td_gaels, 'k*', ...
    EbNoVec, error_td_els, 'ko', ...
    EbNoVec, error_td_omp, 'ro', ...
    EbNoVec, error_td_cosamp, 'rs', ...
    EbNoVec, error_td_mp, 'gs');
legend('LS', 'MMSE', 'GA-MP-LS', 'MP-LS', 'OMP', 'COSAMP', 'MP');
semilogy(EbNoVec, error_td_ls, 'k-', EbNoVec, error_td_mmse, 'b-', EbNoVec, error_td_gaels, 'k-', EbNoVec, error_td_els, 'k-', EbNoVec, error_td_omp, 'r-', EbNoVec, error_td_cosamp, 'r-', EbNoVec, error_td_mp, 'g-');
hold off

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
    fileName = sprintf('MU_MMIMO_M%s_K%s_Np%s_rake_like_detector_mmse_combiner_pilots_%s_v%d_%s.fig',m_rx_antennas,m_tx_antennas,n_pilots,pilot_type,script_version,timeStamp);
    savefig(figura,fileName);
    clear figura
    
    fileName = sprintf('Est_Error_MU_MMIMO_M%s_K%s_Np%s_rake_like_detector_mmse_combiner_pilots_%s_v%d_%s.fig',m_rx_antennas,m_tx_antennas,n_pilots,pilot_type,script_version,timeStamp);
    savefig(fdee_figure,fileName);
    clear fdee_figure
    
    % Save workspace to MAT-file.
    fileName = sprintf('MU_MMIMO_M%s_K%s_Np%s_rake_like_detector_mmse_combiner_pilots_%s_v%d_%s.mat',m_rx_antennas,m_tx_antennas,n_pilots,pilot_type,script_version,timeStamp);
    save(fileName);
end