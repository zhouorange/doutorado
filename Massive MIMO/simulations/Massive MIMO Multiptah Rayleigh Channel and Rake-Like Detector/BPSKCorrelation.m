clear all;close all;clc

K = 1;

NCP = 512;

NFFT = 1e6;

% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*NFFT;                                                                                    % Number of symbols, i.e., number of terminals. The number of pilots must be left out.

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');

Tx = zeros(2,NFFT);
for k=1:1:2
    
    % TX Data.
    % Create array of bits to modulate.
    msg = randi(hStr, [0 1], modOrd, numSym);
    
    % Modulate data.
    source = modulate(hMod, msg);
    
    % Split source among K terminals.
    Tx(k,:) = reshape(source, K, numel(source)/K);
    
end

teste = sum(Tx(1,:).*conj(Tx(2,:)))/NFFT;

a=1;