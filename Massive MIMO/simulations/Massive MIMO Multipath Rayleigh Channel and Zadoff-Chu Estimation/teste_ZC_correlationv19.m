clear all;clc;close all;

show_figures = false;

%% ---------- PRACH Definitions ----------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
NSUBFRAME = 30720;
Ncp = 3168;
NG = 2976;
%v = [0 5 10 15 20 25 30 35 40 45];
v = [0 5 10 15 20 25 30 35 40 45];
Ncs = 13;
position = mod(Nzc-v.*Ncs,Nzc) + 1;                                 % Position of the start of the window. Plus one to correct address matlab vectors.
prach_offset = 10;

%% ---------- Multiptah Channel defintions ----------
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 1000;                                                              % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
Ts = 1/(15000*NFFT);                                                % System Sampleing Rate.

% Small scale fading.
delay = [0 0.977]*1e-6;                                             % Delay in microseconds.
gain  = [0 0];                                                      % Gain in dB (indoor).
numPaths = length(delay);                                           % Number of paths per channel.
totalNumPaths = M*K*numPaths;                                       % Total number of paths between the various antennas and Base Stations.

%% ---------- Initializations ----------
signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);
preambles = complex(zeros(K,NSUBFRAME),zeros(K,NSUBFRAME));

%% ------- Generate Root Zadoff-Chu sequence. -------
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);

%% ****************************** PRACH Transmission ******************************
for i=1:1:K
    
    % ------- Apply Cyclic Shift to Root Zadoff-Chu sequence. -------
    Cv = v(i)*Ncs;
    xuv = xu_root(mod((n+Cv),Nzc)+1);
    
    % ------- Apply DFT to the Preamble. -------
    Xuv = fft(xuv,Nzc);
    
    % ------- Subcarrier Mapping. -------
    bb_signal = [complex(zeros(1,prach_offset),zeros(1,prach_offset)), Xuv, complex(zeros(1,NIDFT-prach_offset-Nzc),zeros(1,NIDFT-prach_offset-Nzc))];
    
    % ------- Apply IDFT to the Baseband Signal. -------
    prach = ifft(bb_signal,NIDFT);
    
    % ------- Add CP. -------
    prach_cp = [prach(NIDFT-Ncp+1:NIDFT), prach];
    
    % ------- Add Guard-time (GT) interval. -------
    y = [prach_cp, zeros(1,NG)];
    
    if(show_figures)
        figure(4);
        stem(0:1:NIDFT-1,abs(fft(y(Ncp+1:NIDFT+Ncp),NIDFT)));
        title('Transmitted PRACH signal.');
    end
    
    preambles(i,:) = y;
      
end

%% *************************** Setup Multipath Massive MIMO Channel ***************************

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

rng(55);
x = [(randn((LS-1)/2+1,totalNumPaths,'double') + 1i*randn((LS-1)/2+1,totalNumPaths,'double')); zeros(N_chann-LS,totalNumPaths); (randn((LS-1)/2,totalNumPaths,'double')) + 1i*randn((LS-1)/2,totalNumPaths,'double')]/sqrt(2);
ch = ifft(x .* repmat(sqrt(S),1,totalNumPaths)) * N_chann / sqrt(LS);

% Plot doppler spectrum figures.
if(show_figures)
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

% ************ Apply power delay profile and large-scale fading to the multiptah channel matrix, i.e., small-scale fading. ************
for idx_ch = 1 : N_chann
    
    % Atualizacao das matrizes de canal para a celula alvo (l = 1): existe uma matriz de canal por percurso do terminal.
    G = reshape(ch(idx_ch,:), M, K, length(pos));
    G(:,:,1) = (g(pos(1)) * G(:,:,1));
    for k = 2:length(pos)
        G(:,:,k) = (g(pos(k)) * G(:,:,k));
    end
    ch(idx_ch,:) = reshape(G, 1, totalNumPaths);
end

%---------- Multipath Channel plus Noise ----------
idx_ch = 1;
H = reshape(ch(idx_ch,:), M, K, length(pos));

y_channel = H(:,:,1)*preambles;
for k = 2:length(pos)
    aux = [complex(zeros(M,(pos(k)-1)),zeros(M,(pos(k)-1))) H(:,:,k)*preambles(:,1:end-(pos(k)-1))];
    y_channel = y_channel + aux;
end

%% ****************************** PRACH Reception ******************************

% ------- CP and GT Removal. -------
rec_signal = y_channel(:,Ncp+1:NIDFT+Ncp);

if(show_figures)
    figure(5);
    signal = fft(rec_signal,NIDFT,2);
    stem(0:1:NIDFT-1,abs(signal(1,:)));
    title('Received base-band signal');
end

% ------- Apply DFT to received signal. -------
rec_fft = fft(rec_signal,NIDFT,2);

% ------- Sub-carrier de-mapping. -------
rec_Xuv = rec_fft(:,prach_offset+1:prach_offset+Nzc);

% ------- Apply DFT to Root Zadoff-Chu sequence. -------
Xu_root = fft(xu_root, Nzc);

% ------- Multiply Local Zadoff-Chu root sequence by received sequence. -------
conj_Xu_root = conj(Xu_root);
multiplied_sequences = complex(zeros(M,Nzc),zeros(M,Nzc));
for mm=1:1:M
    multiplied_sequences(mm,:) = (rec_Xuv(mm,:).*conj_Xu_root);
end

% ------- Squared modulus used for peak detection. -------
NIFFT_CORR = 839;
pdp_freq = ifft(multiplied_sequences,NIFFT_CORR,2)/Nzc;
pdp_freq_adjusted = pdp_freq;

%% --------------------------------------------------------------------------
if(1) % Apenas um pico por vez
    
    PDP_ADJUST_FACTOR = [1, 0.989338652615384 + 0.152095818347174i, 0.957498745916891 + 0.301528270811912i, 0.904911629047376 + 0.445660126688108i, 0.050047412688511 + 1.101839958815881i, -0.473961715008805 + 1.083869828182508i];
    delay_zcorr1 = (pos(numPaths) - 1)/30;
    
    for mm=1:1:M
        for idx=1:1:length(position)         
            pdp_freq_adjusted(mm,position(idx)) = PDP_ADJUST_FACTOR(1)*pdp_freq_adjusted(mm,position(idx));
            if(numPaths > 0)
                pdp_freq_adjusted(mm,position(idx)+delay_zcorr1) = (PDP_ADJUST_FACTOR(2))*pdp_freq_adjusted(mm,position(idx)+delay_zcorr1);
            end
        end
    end
    
    % ----------- Channel estimation. -----------
    H_estimated = zeros(M,K,numPaths);
    for idx=1:1:numPaths
        H_estimated(:,:,idx) = complex(zeros(M,K),zeros(M,K));
    end
    
    for mm=1:1:M
        for kk=1:1:K
            for idx=1:1:numPaths
                H_estimated(mm,kk,idx) = pdp_freq_adjusted(mm,(position(kk)+(idx-1)));
            end
        end
    end
    % -----------------------------------------
    
    error1 = 0;
    for jj=1:1:length(delay)
        for mm=1:1:M
            for kk=1:1:K
                error = abs(H(mm,kk,jj)-pdp_freq_adjusted(mm,position(kk)+((pos(jj)-1)/30)));
                error1 = error1 + error;
                if(error > 0.07)
                    fprintf(1,'h(%d,%d,%d): %s - pdp_freq_adjusted(%d,%d): %s - error: %f\n',mm,kk,jj,num2str(H(mm,kk,jj)),mm,(position(kk)+((pos(jj)-1)/30)),num2str(pdp_freq_adjusted(mm,position(kk)+((pos(jj)-1)/30))),error);         
                end
            end
        end
    end
    error1 = error1 / (M*K*numPaths);
end

%% -------------------------------------------------------------------------
pdp = abs(pdp_freq_adjusted).^2;

if(0)
    for mm=1:1:M
        figure;
        stem(0:1:NIFFT_CORR-1,pdp(mm,:))
        
        figure;
        pdp_freq_adjusted_inv=1./pdp_freq_adjusted(mm,:);
        stem(0:1:NIFFT_CORR-1,abs(pdp_freq_adjusted_inv))
    end
end

error2 = 0;
for idx=1:1:numPaths
    error_matrix = (abs(H - H_estimated));
    error2 = error2 + sum(sum((error_matrix(:,:,idx))));
end
error2 = error2 / (M*K*numPaths);

error1
error2


