clear all;close all;clc

K = 10;

N = 1e5;

NFFT = 2;                                                                                        % Number of points used by the OFDM.
AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*NFFT;

% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');
hDemod = modem.pskdemod(hMod);

Tx = zeros(K,NFFT,N);
for i=1:1:N
    
    % TX Data.
    % Create array of bits to modulate.
    msg = randi(hStr, [0 1], modOrd, numSym);
    
    % Modulate data.
    source = modulate(hMod, msg);
    
    % Split source among K terminals.
    Tx(:,:,i) = reshape(source, K, NFFT);
end

xc = zeros(K,K,NFFT);
coluna = zeros(N,1);
linha = zeros(1,N);
for sc=1:1:NFFT
    for l=1:1:K 
        for c=1:1:K
            
            linha(1,:)  = Tx(l,sc,:);
            coluna(:,1) = Tx(c,sc,:);
            xc(l,c,sc)  = (linha*conj(coluna))/N;
            
        end
    end
end