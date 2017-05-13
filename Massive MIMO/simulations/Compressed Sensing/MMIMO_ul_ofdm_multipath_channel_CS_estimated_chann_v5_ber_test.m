clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 32;                                                                                            % Number of pilots per OFDM symbol.

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
EbNoVec = -20:2:20;                                                                                 % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10((NFFT-K*Np+Np)/NFFT) + 10*log10(modOrd);    % converting to symbol to noise ratio.
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                                                         % Calculate SNR from EsNo in dB.

nTotalOfBits = 1e8;
nErrors = 1000000;
                                                                                  
%% @@@@@@@@@@@@@@@@ Create Pilots @@@@@@@@@@@@@@@@@@@@@@@@@@@@
fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:Np-1)/Np);              % Normalized pilot signal, i.e., unit power.

Pdiag = diag(P);

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
[BER_MRC, BER_ZF_LE, BER_MMSE_LE, BER_MFB, BER] = deal(zeros(1, length(EbNoVec)));

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% Set up a figure for visualizing BER results.
h = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-7 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
title('Massive MU-MIMO on Uplink');

%------------ Create Pilots ------------
if(Np==16) % Funciona somente com OMP.
    [ppos, flag_data] = getOptimumPpos16();
elseif(Np==32) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos32();
elseif(Np==34) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos34();
elseif(Np==35) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos35();
elseif(Np==36) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos36();
elseif(Np==40) % Valores ótimos para LS e OMP.
    [ppos, flag_data] = getOptimumPpos40();
else
    error('Invalid value for Np!!!!!');
end

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
for idx = 1:length(snr)
    
    linearSnr = 10^(0.1*snr(idx));
    
    ofdm_symbol_number = 0;
    idx_ch = 1;
    x = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    aux = complex(zeros(K,NFFT+NCP),zeros(K,NFFT+NCP));
    
    nErrs = 0;
    nBits = 0;
    iter = 0;
    
    while((nErrs < nErrors) && (nBits < nTotalOfBits))
        
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
            for c_idx=1:1:Np
                Tx(l_idx,ppos(l_idx,:)) = P;
            end
        end
        
        % Create OFDM symbol.
        sequence = (NFFT/sqrt(NFFT-Np*K+Np))*ifft(Tx,NFFT,2);
        
        % Add CP.
        ofdm = [sequence(:,NFFT-NCP+1:end), sequence];
        
        %--------------------------- Channel ------------------------------
        if(add_awgn)
            % Add channel noise power to faded data.
            r = awgn(ofdm, snr(idx), 0, hStr);
        else
            r = ofdm;
        end
        
        %------------------- Reception (base station) ---------------------
        
        % Remove CP.
        rx = r(:,NCP+1:end);
        
        % OFDM.
        ofdm_rx = (sqrt(NFFT-Np*K+Np)/NFFT)*fft(rx,NFFT,2);
        
        % Retrieve modulation symbols.
        data_symbols = ofdm_rx(~flag_data).';
        
        % Demodulate signal.
        demoded = demodulate(hDemod, data_symbols);
        
        %% ----------- Collecting Errors ---------------
        nErrs = nErrs + biterr(msg, demoded);
        
        nBits = nBits + length(msg(:));
        
        fprintf(1,'SNR: %d\n',EbNoVec(idx));
        fprintf(1,'BER: %f - nErrs: %d - nBits: %d - iter: %d\n',(nErrs./nBits),nErrs,nBits,iter);
        
    end
    
    % Calculate BER for current point
    BER(idx) = nErrs./nBits;
    
end

% Plot results
semilogy(EbNoVec(1:idx), BER(1:idx), 'ro');
berTheory = berawgn(EbNoVec(1:idx),'psk',AlphabetSize,'nondiff');
semilogy(EbNoVec,berTheory,'k')
legend('Simulation','Theory')
hold off;
