clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 37;                                                                                            % Number of pilots per OFDM symbol. (use prime numbers)

L = 1;                                                                                              % Number of cells.
K = 1;                                                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 1;                                                                                            % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*(NFFT - K*Np);                                                                           % Number of symbols, i.e., number of terminals. The number of pilots must be leaved out.
NCP = 512;                                                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
delta_f = 15000;                                                                                    % Frequency space between subcarriers.
Ts = 1/(delta_f*NFFT);                                                                              % System Sampleing Rate.
numSymbInSubframe = 12;                                                                             % Number of symbols in a subframe (1ms). 12 OFDM symbols for extended CP.

add_awgn = false;                                                                                    % Enable/disable addtion of additive gaussian noise to signals.
EbNoVec = 1500; %-4:2:6; %-5:1:5;%-20:2:40;                                                                                  % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10((NFFT-K*Np+Np)/NFFT) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                                                  % Radius given in meters.
cellHole = 100;                                                                                     % Cell hole in meters.

power_ctrl_enabled = true;                                                                          % enable/disable power control at eNodeB.

nTotalOfTrials = 1e5;
nTotalOfBits = 1e7;
nErrors = 1e6;
debug = false;

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 30]*Ts;                                             % Delay in microseconds.
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
    fprintf(1,'Invalid value for Np!!!!!');
end

if(K==1)
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
[BER_MMSE_ELS] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
if(plot_fig)
    figura = figure; grid on; hold on;
    set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
    xlabel('Eb/No (dB)'); ylabel('BER'); set(figura,'NumberTitle','off');
    set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
    strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
    title(strTitle);
end

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
a = (NFFT/sqrt(NFFT-(Np*K)+Np));
N_chann = 10;
ChannelOrder = pos(length(pos))-1;
mp_correct_estimation = zeros(1,length(snr));
omp_correct_estimation = zeros(1,length(snr));
for snr_idx = 1:1:length(snr)
    
    linearSnr = 10^(snr(snr_idx)/10);
    
    ofdm_symbol_number = 0;
    idx_ch = 1;
    iter = 0;
    nTrial = 0;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    nErrs_mmse_els = 0;
    
    while(nTrial < nTotalOfTrials)
        
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
        
        x = [complex(zeros(M,(pos(1)-1)),zeros(M,(pos(1)-1))) H(:,:,1)*ofdm complex(zeros(M,(pos(length(pos))-pos(1))),zeros(M,(pos(length(pos))-pos(1))))];
        for k = 2:length(pos)
            aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*ofdm complex(zeros(M,(pos(length(pos))-pos(k))),zeros(M,(pos(length(pos))-pos(k))))];
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
            H_hat_els = zeros(M, K, length(pos));
            for m_idx=1:1:M
                for k_idx=1:1:K
                    
                    Pdiag = diag(P(k_idx,:));
                    
                    %% ********* Fourier Basis *********
                    Fl = F(ppos(k_idx,:).',:);
                    
                    y = ofdm_rx(m_idx,ppos(k_idx,:)).';
                    
                    
                    %-------------------
                    %P(k_idx,:) = exp((-1i.*pi.*839.*n.*(n+1))./Np);
%                     test = ifft(y.*conj(P(k_idx,:).'),Np);
%                     
%                     stem(abs(test));
                    
                    %sum_test = sum(abs(test));
                    
                    % --------- Iterative Greedy Algorithm ---------
                    A = Pdiag*Fl;
                    
%                     % OMP.
%                     g_hat_omp = OMP_origv1(A,y,numPaths);
%                     
%                     if((g_hat_omp(1)~=0 && g_hat_omp(31)~=0))
%                         omp_correct_estimation(snr_idx) = omp_correct_estimation(snr_idx) + 1;
%                     end
                    
                    
                    
                    % MP.
                    [g_hat_mp tap_pos tap_res]= mp(A,y);
                    
%                     if((tap_pos(1)==1 && tap_pos(2)==31) || (tap_pos(1)==31 && tap_pos(2)==1))
%                         mp_correct_estimation(snr_idx) = mp_correct_estimation(snr_idx) + 1;
%                     end

%                     figure
%                     stem(abs(g_hat_mp))
%                     figure
%                     stem(tap_res)
%                     
%                     aaaa=1;
                    
                    % --------- Linear Algorithms ---------
                    
                    % E-LS: Enhanced Least Squares.
                    % Forma que melhora a performance do MMSE, porém, deve-se saber a posição dos taps diferentes de zero.
                    
                    
                    %d = (((Fl(:,1:pos(length(pos)))')*(Pdiag')*(Pdiag*Fl(:,1:pos(length(pos)))))^(-1))*(Fl(:,1:pos(length(pos)))' * Pdiag')*y;
                    
                    %d_abs = abs(d);
                    
                    %fprintf(1,'min_tap_pos(1): %d - min_tap_pos(2): %d - min_error: %d \n',min_tap_pos(1),min_tap_pos(2),min_error);
                    
                    
                    %
                    %                     d = zeros(1,pos(length(pos)));
                    %                     d(min_tap_pos) = [1 1];
                    %                     D = diag(d);
                    %                     B = Pdiag*Fl(:,1:pos(length(pos)))*D;
                    %
                    %
                    %                     Bl = B(:,[min_tap_pos(1) min_tap_pos(length(min_tap_pos))]);
                    %                     invMatEls = (((Bl'*Bl)^(-1))*Bl');
                    %                     g_hat_els = invMatEls*y;
                    %                     g_hat_els = [g_hat_els(1); zeros(31-2,1); g_hat_els(2)];
                    %
                    %                     % --------- Estimate of H ---------
                    %                     H_hat_els(m_idx,k_idx,:) = g_hat_els(min_tap_pos);
                    
                    nTrial = nTrial + 1;
                end
            end
        end
        %------------------------------------------------------------------
        
        % Slicer.
        yol = complex(zeros(1,31),zeros(1,31));
        summ = ofdm_rx(1);
        ofdm_rx_ext = [ofdm_rx complex(zeros(1,31),zeros(1,31))];
        jj = 31;
        yo = 0;
        for ii=2:1:length(ofdm_rx_ext)
            
            yol(jj) = yol(jj-1);
            yol(1) = yo;
            
            jj = jj + 1;
            if(jj > 31)
                jj = 2;
            end
            
            yob = demodulate(hDemod, summ);
            yo = modulate(hMod, yob);
            
            
            
            
            summ = ofdm_rx_ext(ii) - H(1,1,2)*yol(31);
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
    fprintf(1,'SNR: %d - mp_correct_estimation: %d - nTrial: %d - ratio: %d \n',EbNoVec(snr_idx),mp_correct_estimation(snr_idx),nTrial,(mp_correct_estimation(snr_idx)/nTrial));
    fprintf(1,'SNR: %d - omp_correct_estimation: %d - nTrial: %d - ratio: %d\n',EbNoVec(snr_idx),omp_correct_estimation(snr_idx),nTrial,(omp_correct_estimation(snr_idx)/nTrial));
end

figure
hold on
plot(EbNoVec,mp_correct_estimation)
plot(EbNoVec,omp_correct_estimation)
xlabel('Eb/No (dB)'); ylabel('Correct Estimation'); 
grid on
hold off

mp_correct_estimation_ratio = (mp_correct_estimation./nTrial);
omp_correct_estimation_ratio = (omp_correct_estimation./nTrial);
figure
hold on
plot(EbNoVec,mp_correct_estimation_ratio)
plot(EbNoVec,omp_correct_estimation_ratio)
xlabel('Eb/No (dB)'); ylabel('Correct Estimation Ratio'); 
grid on
hold off


