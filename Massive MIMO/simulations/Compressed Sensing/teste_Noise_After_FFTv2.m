clear all;close all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                                                   % Disable/enable figure plotting.

Np = 101;                                                                                           % Number of pilots per OFDM symbol. (use prime numbers)

L = 1;                                                                                              % Number of cells.
K = 1;                                                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
NFFT = 2048;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*(NFFT - K*Np);                                                                           % Number of symbols, i.e., number of terminals. The number of pilots must be leaved out.
NCP = 512;                                                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.


% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

%% ----------------------------------- Set up the simulation -----------------------------------
% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

%% @@@@@@@@@@@@@@@@ Create Pilots @@@@@@@@@@@@@@@@@@@@@@@@@@@@
fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:Np-1)/Np);              % Normalized pilot signal, i.e., unit power.

Pdiag = diag(P);

F = fft(eye(NFFT));
F = F(:,1:NCP);

%------------ Retrieve Pilot Positions ------------
ppos = randperm(NFFT,Np*K);
ppos = reshape(ppos,K,Np);
ppos = sort(ppos,2);

flag_data = false(K,NFFT);
for l_idx=1:1:K
    for c_idx=1:1:Np
        flag_data(:,ppos(l_idx,c_idx)) = true(K,1);
    end
end

%% ----------------------------------- Loop over selected EbNo points. -----------------------------------
NumIter = 10000;
var_sequence = 0;
var_ofdm_rx = 0;
for iter_idx = 1:1:NumIter
    
    %---------- Transmission (UE) ----------
    % Create array of bits to modulate.
    msg = randi(hStr, [0 1], modOrd, numSym);
    
    % Modulate data.
    source = modulate(hMod, msg);
    
    % Split source among K terminals.
    Tx = complex(zeros(K,NFFT),zeros(K,NFFT));
    tt = flag_data(1:K,:);
    Tx(~tt) = reshape(source, K, numel(source)/K); clear source;
    for l_idx=1:1:K
        Tx(l_idx,ppos(l_idx,:)) = P(l_idx,:);
    end
    
    var_p = var(P);
    
    var_tx = var(Tx);
    
    % Create OFDM symbol.
    sequence = (NFFT/sqrt(NFFT-(Np*K)+Np))*ifft(Tx,NFFT,2);
    
    for ii=1:1:K
        var_sequence = var_sequence + var(sequence(ii,:));
    end
    
    % Add CP.
    ofdm = [sequence(:,NFFT-NCP+1:end), sequence];
    
    %------------------- Reception (base station) ---------------------
    rx = ofdm;
    
    % Remove CP.
    rx = rx(:,NCP+1:end);
    
    % Retrieve modulation symbols by applying FFT to received signal.
    ofdm_rx = (sqrt(NFFT-(Np*K)+Np)/NFFT)*fft(rx,NFFT,2);
    
    for ii=1:1:K
        var_ofdm_rx = var_ofdm_rx + var(ofdm_rx(ii,:));
    end
    
end

var_ofdm_rx = var_ofdm_rx/NumIter;
var_sequence = var_sequence/NumIter;

a=1;