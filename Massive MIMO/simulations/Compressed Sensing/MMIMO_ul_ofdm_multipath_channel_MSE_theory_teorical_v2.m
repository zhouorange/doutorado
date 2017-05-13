clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 37;                                                                                           % Number of pilots per OFDM symbol. (use prime numbers)

L = 1;                                                                                              % Number of cells.
K = 1;                                                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 2;                                                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*(NFFT - K*Np);                                                                           % Number of symbols, i.e., number of terminals. The number of pilots must be leaved out.
NCP = 512;                                                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = true;                                                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 100:2:102;                                                                                % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10((NFFT-K*Np+Np)/NFFT) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                                                  % Radius given in meters.
cellHole = 100;                                                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e8;
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
Fs_chann = 500;                                                     % Channel Sampling Rate in Hz. The channel is sampled every 2 ms. (Periodo de amostragem do canal ~1 amostra / Slot OFDM)
Ts_chann = 1/Fs_chann;
NumIter = 256;                                                      % Number of samples used to sample the channel. Duration of the channel in  number of samples.
delta_f_chann = Fs_chann/NumIter;                                   % in Hz.
f = -Fs_chann/2:delta_f_chann:Fs_chann/2;
snr_idx = find(f<fd);
f = f(snr_idx);
snr_idx = find(f>-fd);
f = f(snr_idx);
f = f.';

LS = length(f);
S = 1/pi/fd./sqrt(1 - (f/fd).^2) * 10^(Pd/10);
S = S * LS / sum(S);                                                % Energy nomalization.
S1 = S;
S = [S((LS+1)/2:LS); zeros(NumIter-LS,1); S(1:(LS-1)/2)];           % Interpolation of the Doppler Spectrum by N_chann/LS.

% ************************* Generate Multipath Massive MIMO Channel. **********************
rng(55);
x = [(randn((LS-1)/2+1,totalNumPaths,'double') + 1i*randn((LS-1)/2+1,totalNumPaths,'double')); zeros(NumIter-LS,totalNumPaths); (randn((LS-1)/2,totalNumPaths,'double')) + 1i*randn((LS-1)/2,totalNumPaths,'double')]/sqrt(2);
ch = ifft(x .* repmat(sqrt(S),1,totalNumPaths)) * NumIter / sqrt(LS);

% Plot doppler spectrum figures.
if(plot_fig)
    figure;
    plot(f, abs(S1));
    
    figure;
    plot(S)
    
    figure;
    Tsd = 1/Fs_chann;
    TNd = Tsd*NumIter; % Channel duration in seconds.
    plot( (0:NumIter-1)*Tsd, 10*log10(abs(ch)) );
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
if(plot_fig)
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
end

% Calculate path-loss for all users, in meters.
path_loss = radius.^gamma;

% ************ Apply power delay profile and large-scale fading to the multiptah channel matrix, i.e., small-scale fading. ************
for iter = 1 : NumIter
    
    if(power_ctrl_enabled)
        largeScaleFading = 1;
    else
        % Calculate shadowing for each one of different channels.
        shadowing_att = lognrnd(0,sshadow,1,K);
        % Large scale fading calculated according to Marzetta.
        largeScaleFading = repmat(sqrt(shadowing_att./path_loss), M, 1);
    end
    
    % Atualizacao das matrizes de canal para a celula alvo (l = 1): existe uma matriz de canal por percurso do terminal.
    G = reshape(ch(iter,:), M, K, length(pos));
    G(:,:,1) = (g(pos(1)) * G(:,:,1)) .* largeScaleFading;
    for k = 2:length(pos)
        G(:,:,k) = (g(pos(k)) * G(:,:,k)) .* largeScaleFading;
    end
    ch(iter,:) = reshape(G, 1, totalNumPaths);
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

IF = (sqrt(NFFT))*ifft(eye(NFFT));

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

if(K~=10)
    flag_data = [false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false true false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false true false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false true true false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false true true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false true false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false false];
    ppos = [54 95 116 141 218 225 231 267 363 377 463 505 530 560 617 712 724 890 968 983 998 1035 1065 1197 1226 1227 1249 1291 1300 1301 1487 1502 1649 1680 1749 1827 1899];
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
[BER_MMSE_LE, BER_MFB] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-8 1000]);
xlabel('Eb/No (dB)'); ylabel('MSE'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
title(strTitle);

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
H = reshape(ch(1,:), M, K, length(pos));

a = (NFFT/sqrt(NFFT-(Np*K)+Np));

mmse_error_theory = zeros(M,K,length(snr));
ls_error_theory = zeros(M,K,length(snr));
mse_siso_ls = zeros(M,K,length(snr));
mse_siso_mmse = zeros(M,K,length(snr));
error_omp = zeros(M,K,length(snr));
error_ls = zeros(M,K,length(snr));
error_mmse = zeros(M,K,length(snr));
NumIter = 1000;
ChannelOrder = pos(length(pos))-1;
for snr_idx = 1:1:length(snr)
    
    linearSnr = 10^(snr(snr_idx)/10);
    
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
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    x = [H(:,:,1)*ofdm complex(zeros(M,(pos(length(pos))-1)),zeros(M,(pos(length(pos))-1)))];
    for k = 2:length(pos)
        aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*ofdm];
        x = x + aux;
    end
    
    for iter=1:1:NumIter
        
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
        
        % Massive MIMO Channel Estimation.
        h = zeros(pos(length(pos)),1);
        H_hat_omp = zeros(M, K, length(pos));
        H_hat_ls = zeros(M, K, length(pos));
        H_hat_mmse = zeros(M, K, length(pos));
        for m_idx=1:1:M
            for k_idx=1:1:K
                
                Pdiag = diag(P(k_idx,:));
                
                %% ********* Fourier Basis *********
                Fl = F(ppos(k_idx,:).',:);
                
                IFl = IF(1:37,ppos(k_idx,:));
                
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
%                 invMatMmse = ((A*A' + ((NFFT*numPaths)/(linearSnr*(a.^2)))*eye(Np)) )^-1;
%                 invMatMmse = A'*invMatMmse;
%                 g_hat_mmse = invMatMmse*y;

%                 % Implementação através da versão no dominio da frequência de h.
%                 invMatMmse = ((Pdiag*Pdiag' + ((NFFT)/(linearSnr*(a.^2)))*eye(Np)) )^-1;
%                 invMatMmse = Pdiag'*invMatMmse;
%                 G_hat_mmse = invMatMmse*y;                
%                 g_hat_mmse = ifft(G_hat_mmse,31);
%                 g_hat_mmse = g_hat_mmse(1:pos(length(pos)));
% 
                % Implementação através a versão no dominio da frequência de h.
                G_hat_mmse = complex(zeros(NFFT,1),zeros(NFFT,1));
                invMatMmse = ((Pdiag'*Pdiag + ((NFFT)/(linearSnr*(a.^2)))*eye(Np)) )^-1;
                invMatMmse = invMatMmse*Pdiag';
                G_hat_mmse(ppos(k_idx,:),:) = invMatMmse*y;
                g_hat_mmse = (2048/37)*ifft(G_hat_mmse,NFFT);
                g_hat_mmse = g_hat_mmse(1:pos(length(pos)));

                
                % Interpolar G_hat_mmse com a função: interp1 antes de
                % aplicar a IFFT.                
                interp_data = invMatMmse*y;
                interp_data = interp_data.';
                interp_pos = ppos(k_idx,:);                              
                vq = interp1(interp_pos, interp_data, 1:1:2048,'spline','extrap');
                
                H_interpolated = interpolate(interp_data,interp_pos,NFFT,'spline');                
                g_hat_mmse_interp = ifft(H_interpolated,NFFT);
                g_hat_mmse_interp = g_hat_mmse_interp(1:pos(length(pos)));

%                 % Forma que melhora a performance do MMSE, porém, deve-se saber a posição dos taps diferentes de zero.
%                 d = zeros(1,pos(length(pos)));
%                 d(pos) = [1 1];
%                 D = diag(d);
%                 B = Pdiag*Fl(:,1:pos(length(pos)))*D;                 
%                 invMatMmse = (B*B' + ((NFFT*numPaths)/(linearSnr*(a.^2)))*eye(Np))^-1;
%                 invMatMmse = B'*invMatMmse;
%                 g_hat_mmse = invMatMmse*y;
                
                % --------- Estimate of H ---------
                h(pos,1) = [H(m_idx,k_idx,1);H(m_idx,k_idx,2)];
                
                h_nfft = [h.' complex(zeros(1,NFFT-length(h)),zeros(1,NFFT-length(h)))];
                H_fft = fft(h_nfft);
                
                error_omp(m_idx,k_idx,snr_idx) = error_omp(m_idx,k_idx,snr_idx) + (sum(abs(h - g_hat_omp(1:pos(length(pos)))).^2)/pos(length(pos)));
                error_ls(m_idx,k_idx,snr_idx) = error_ls(m_idx,k_idx,snr_idx) + (sum(abs(h - g_hat_ls).^2)/pos(length(pos)));
                error_mmse(m_idx,k_idx,snr_idx) = error_mmse(m_idx,k_idx,snr_idx) + (sum(abs(h - g_hat_mmse).^2)/pos(length(pos)));

                
                
                
                n = awgn(complex(zeros(NFFT,1),zeros(NFFT,1)), snr(snr_idx), 0, hStr);
                N = (sqrt(NFFT-(Np*K)+Np)/NFFT)*fft(n,NFFT);
                N = N(ppos(k_idx,:),1);
                
                % LS Error Theory.
                error = ((Fl(:,1:pos(length(pos)))'*Fl(:,1:pos(length(pos))))^-1)*(Fl(:,1:pos(length(pos)))'*Pdiag')*N;
                error = (error'*error)/pos(length(pos));
                ls_error_theory(m_idx,k_idx,snr_idx) = ls_error_theory(m_idx,k_idx,snr_idx) + error;

                % MMSE Error Theory
                b = (NFFT*numPaths)/((a.^2)*linearSnr);              
                InvMatrix = (A*A' + b*eye(Np))^-1;
                h_hat_mmse = A'*InvMatrix*(A*h + N);
                error = h - h_hat_mmse;
                error = (error'*error)/pos(length(pos));
                mmse_error_theory(m_idx,k_idx,snr_idx) = mmse_error_theory(m_idx,k_idx,snr_idx) + error;


                if(iter==1)
                    
                    mse_siso_mmse(m_idx,k_idx,snr_idx) = (NFFT/(linearSnr*(a.^2)))*(trace( ( ( A*A' + ((NFFT*numPaths)/(linearSnr*(a.^2)))*eye(Np) ) )^-1  ) / Np   );

                    mse_siso_ls(m_idx,k_idx,snr_idx) = (NFFT/(linearSnr*(a.^2)))*(trace(((Fl(:,1:pos(length(pos)))'*Fl(:,1:pos(length(pos))))^-1))/pos(length(pos)));
                end
            end
        end
    end
end

error_omp = error_omp/NumIter;
error_ls = error_ls/NumIter;
error_mmse = error_mmse/NumIter;
ls_error_theory = ls_error_theory/NumIter;
mmse_error_theory = mmse_error_theory/NumIter;

for m_idx=1:1:M
    for k_idx=1:1:K
        
        error_ls_vec(1,:) = error_ls(m_idx,k_idx,:);
        error_mmse_vec(1,:) = error_mmse(m_idx,k_idx,:);
        error_omp_vec(1,:) = error_omp(m_idx,k_idx,:);
        mse_siso_mmse_vec(1,:) = mse_siso_mmse(m_idx,k_idx,:);
        mse_siso_ls_vec(1,:) = mse_siso_ls(m_idx,k_idx,:);
        mmse_error_theory_vec(1,:) = mmse_error_theory(m_idx,k_idx,:);
        
        semilogy(EbNoVec, error_omp_vec, 'ro--');
        semilogy(EbNoVec, error_ls_vec, 'cs--');
        semilogy(EbNoVec, error_mmse_vec, 'k*--');
        semilogy(EbNoVec, mse_siso_ls_vec, 'c');
        semilogy(EbNoVec, mse_siso_mmse_vec, 'k');
        
        semilogy(EbNoVec, mmse_error_theory_vec, 'ks--');
    end
end

legend('OMP','LS','MMSE','MSE LS SISO','MSE MMSE SISO','MMSE Theory');
hold off
