clear all;clc;close all;

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
NSUBFRAME = 30720;
Ncp = 3168;
NG = 2976;
v = [10];
Ncs = 13;

K = 1; % Number of Tranmitters, i.e., UEs or number of single-antenna terminals.

prach_offset = 10;

signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);

preambles = complex(zeros(K,NSUBFRAME),zeros(K,NSUBFRAME));

show_figures = false;

% ------- Generate Root Zadoff-Chu sequence. -------
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

%% *************************** Multipath Rayleigh Channel ***************************
delay_zcorr = 10;
delay = delay_zcorr*30;

h110 = 0.5 + 1i*0.8;
h111 = 0.2 + 1i*0.2;

energy = abs(h110).^2 + abs(h111).^2;
normalization_factor = 1/sqrt(energy);

h110 = h110*normalization_factor;
h111 = h111*normalization_factor;

energy = abs(h110).^2 + abs(h111).^2;

if(delay > 0)
    preambles_delayed = [complex(zeros(1,delay),zeros(1,delay)), preambles(1,1:end-delay)];    
    y_channel = (h110*preambles + h111*preambles_delayed);
else
    y_channel = h110*preambles;
end

% preambles_delayed = [complex(zeros(1,delay),zeros(1,delay)), preambles(1,1:end-delay)]; 
% y_channel = preambles_delayed;%h110*preambles_delayed;

%% ****************************** PRACH Reception ******************************

% ------- CP and GT Removal. -------
rec_signal = y_channel(Ncp+1:NIDFT+Ncp);

if(show_figures)
    figure(5);
    stem(0:1:NIDFT-1,abs(fft(rec_signal,NIDFT)));
    title('Received base-band signal');
end

% ------- Apply DFT to received signal. -------
rec_fft = fft(rec_signal,NIDFT);

% ------- Sub-carrier de-mapping. -------
rec_Xuv = rec_fft(prach_offset+1:prach_offset+Nzc);

% ------- Apply DFT to Root Zadoff-Chu sequence. -------
Xu_root = fft(xu_root, Nzc);

% ------- Multiply Local Zadoff-Chu root sequence by received sequence. -------
conj_Xu_root = conj(Xu_root);
multiplied_sequences = (rec_Xuv.*conj_Xu_root);

% ------- Squared modulus used for peak detection. -------
PDP_ADJUST_FACTOR = 1; % delay 0

%PDP_ADJUST_FACTOR = 0.989338652615384 + 0.152095818347174i; % delay 1
%PDP_ADJUST_FACTOR = 0.989338652615384 + 0.152095818347175i; % delay 1

%PDP_ADJUST_FACTOR = 0.957498745916891 + 0.301528270811912i; % delay 2
%PDP_ADJUST_FACTOR = 0.957498745916891 + 0.301528270811912i; % delay 2

%PDP_ADJUST_FACTOR = 0.904911629047376 + 0.445660126688108i; % delay 3

%PDP_ADJUST_FACTOR = 0.050047412688511 + 1.101839958815881i; % delay 10
%PDP_ADJUST_FACTOR = 0.050047412688511 + 1.101839958815881i; % delay 10

%PDP_ADJUST_FACTOR = -0.473961715008805 + 1.083869828182508i; % delay 13
%PDP_ADJUST_FACTOR = -0.473961715008805 + 1.083869828182508i; % delay 13

%pdp_freq = ifft(multiplied_sequences,Nzc)/Nzc;

NIFFT_CORR = 839;

pdp_freq = ifft(multiplied_sequences,NIFFT_CORR)/Nzc;
pdp_freq_adjusted = PDP_ADJUST_FACTOR*pdp_freq;

pdp_freq_adjusted(710) = pdp_freq_adjusted(710);
if(delay > 0)
    pdp_freq_adjusted(710+delay_zcorr) = (0.050047412688511 + 1.101839958815881i)*pdp_freq_adjusted(710+delay_zcorr);
end

pdp = abs(pdp_freq_adjusted).^2;

stem(0:1:NIFFT_CORR-1,pdp)
%stem(0:1:NIFFT_CORR-1,pdp_freq_adjusted)

h110
h111
pdp_freq_adjusted(710)
pdp_freq_adjusted(710+delay_zcorr)


