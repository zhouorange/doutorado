clear all;close all;clc

K = 10;

NFFT = 1e5;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*NFFT; 

% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

% TX Data.
% Create array of bits to modulate.
msg = randi(hStr, [0 1], modOrd, numSym);

% Modulate data.
source = modulate(hMod, msg);

% Split source among K terminals.
Tx = reshape(source, K, numel(source)/K);


xc = (Tx*Tx')/NFFT;

a=1;