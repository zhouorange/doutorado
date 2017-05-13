clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                   % Disable/enable figure plotting.

Np_vector = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151];
Np = 37;                                                            % Number of pilots per OFDM symbol. (Must be a prime number.)

L = 1;                                                              % Number of cells.
K = 1;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 2;                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
modOrd = 2;                                                         % Constellation size = 2^modOrd.
numSym = K*(NFFT - K*Np);                                             % Number of symbols, i.e., number of terminals. The number of pilots must be leaved out.
NCP = 512;                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = false;                                                   % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 10000;                                                    % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                          % enable/disable power control at eNodeB.

nTotalOfBits = 1e8;
nErrors = 1000000;
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

%% @@@@@@@@@@@@@@@@ Create Pilots @@@@@@@@@@@@@@@@@@@@@@@@@@@@
u = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571];
n = (0:1:Np-1);
P = complex(zeros(K,Np),zeros(K,Np));
for u_idx=1:1:K
    P(u_idx,:) = exp((-1i.*pi.*u(u_idx).*n.*(n+1))./Np);              % Normalized pilot signal, i.e., unit power.
end

F = fft(eye(NFFT));
F = F(:,1:NCP);

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

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for idx = 1:length(snr)
    
    linearSnr = 10^(0.1*snr(idx));
    
    ofdm_symbol_number = 0;
    idx_ch = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    flag_data_min = false(K,NFFT);
    ppos_min = zeros(1,M);
    min_error_omp = Inf;
    avg_error_omp = Inf;
    min_error_ls = Inf;
    avg_error_ls = Inf;
    min_error = Inf;
    while(min_error>1e-5)
        
        %------------ Create Pilots ------------
        flag_data = false(K,NFFT);
        ppos = 0;
        while(length(unique(ppos)) ~= Np*K)
            ppos = randperm(NFFT,Np*K);
        end
        ppos = reshape(ppos,K,Np);
        ppos = sort(ppos,2);
        
        for l_idx=1:1:K
            for c_idx=1:1:Np
                flag_data(:,ppos(l_idx,c_idx)) = true(K,1);
            end
        end
        
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
        sequence = (NFFT/sqrt(NFFT-Np*K))*ifft(Tx,NFFT,2);
        
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
            r = awgn(x, snr(idx), 0, hStr);
        else
            r = x;
        end
        
        %------------------- Reception (base station) ---------------------
        % Get a whole subframe.
        rx = r(:,1:end-(pos(length(pos))-1));
        
        % Remove CP.
        rx = rx(:,NCP+1:end);
        
        % Retrieve modulation symbols.
        ofdm_rx = (sqrt(NFFT-Np*K)/NFFT)*fft(rx,NFFT,2);
        
        avg_error_ls = 0;
        avg_error_omp = 0;
        h = zeros(NCP,1);
        for m_idx=1:1:M
            
            for k_idx=1:1:K
                
                Pdiag = diag(P(k_idx,:));
                
                %% ********* Fourier Basis *********
                Fl = F(ppos(k_idx,:).',:);
                
                y = ofdm_rx(m_idx,ppos(k_idx,:)).';
                A = Pdiag*Fl;
                g_hat_omp = OMP_orig(A,y,numPaths);
                
                % H(m,k,n)
                h(pos,1) = [H(m_idx,k_idx,1);H(m_idx,k_idx,2)];
                
                % Error OMP.
                error_omp = norm(h(1:length(g_hat_omp))-g_hat_omp)^2/norm(h(1:length(g_hat_omp)))^2;
                
                % Average error OMP.
                avg_error_omp = avg_error_omp + error_omp;
                
                % ---------Linear Algorithms ---------
                % Least Squares
                A = Pdiag*Fl(:,1:pos(length(pos)));
                invMat = (((A'*A)^(-1))*A');
                g_hat_ls = invMat*y;
                
                % Error LS.
                error_ls = norm(h(1:length(g_hat_ls))-g_hat_ls)^2/norm(h(1:length(g_hat_ls)))^2;
                
                % Average error LS.
                avg_error_ls = avg_error_ls + error_ls;
                
                %                 if(error_cse>1e-10)
                %
                %                     error('too large error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                %                 end
                
            end
        end
        avg_error_omp = avg_error_omp/M*K;
        
        avg_error_ls = avg_error_ls/M*K;
        
        if(Np < 30)
            avg_error = avg_error_omp;
        else
            avg_error = (avg_error_omp + avg_error_ls)/2;
        end
               
        if(avg_error < min_error)
            min_error = avg_error;
            ppos_min = ppos;
            flag_data_min = flag_data;
        end
        
        fprintf(1,'avg error: %d - min error: %d\n',avg_error,min_error);
        
    end
    
end
