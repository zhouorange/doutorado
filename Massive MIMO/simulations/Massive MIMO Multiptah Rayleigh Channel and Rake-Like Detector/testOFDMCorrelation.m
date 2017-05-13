clear all;close all;clc

K = 1;

NCP = 512;

NFFT = 2048;

% Create a local random stream to be used by random number generators for repeatability.
hStr = RandStream('mt19937ar');

AlphabetSize = 4;                                                                                   % Modulation alphabet size (QPSK)
modOrd = log2(AlphabetSize);                                                                        % Bits/symbol
numSym = K*NFFT;                                                                                    % Number of symbols, i.e., number of terminals. The number of pilots must be left out.

% Create QPSK mod-demod objects.
hMod = modem.pskmod('M', 2^modOrd, 'SymbolOrder', 'gray', 'InputType', 'bit');

numIter = 100000;
ofdm = zeros(numIter,NFFT+NCP);
for ii=1:1:numIter
    
    % TX Data.
    % Create array of bits to modulate.
    msg = randi(hStr, [0 1], modOrd, numSym);
    
    % Modulate data.
    source = modulate(hMod, msg);
    
    % Split source among K terminals.
    Tx = reshape(source, K, numel(source)/K);
    
%     varTx = var(Tx);
%     varTx2 = sum(Tx.*conj(Tx))/length(Tx);
%     meanTx = mean(imag(Tx));
    
    % Create OFDM symbol.
    sequence = (NFFT/sqrt(NFFT))*ifft(Tx,NFFT,2);
    
    % Add CP.
    ofdm(ii,:) = [sequence(:,NFFT-NCP+1:end), sequence];
    
end

for kk=1:1:(NFFT+NCP)
    
    
    teste = sum(ofdm(:,kk).*conj(ofdm(:,kk)))/numIter;
    
    a=1;
    
end