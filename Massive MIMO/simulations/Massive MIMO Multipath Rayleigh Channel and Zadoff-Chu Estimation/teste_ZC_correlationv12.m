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
delay_zcorr1 = 2;
delay1 = delay_zcorr1*12;

delay_zcorr2 = 2;
delay2 = delay_zcorr2*12;

h110 = 0.5 + 1i*0.8;
h111 = 0.2 + 1i*0.2;
h112 = 0.7 + 1i*0.3;

energy = (abs(h110).^2 + abs(h111).^2 + abs(h112).^2) / 3;
normalization_factor = 1/sqrt(energy);

h110 = h110*normalization_factor;
h111 = h111*normalization_factor;
h112 = h112*normalization_factor;

energy = (abs(h110).^2 + abs(h111).^2 + abs(h112).^2) / 3;

if(delay1 > 0 || delay2 > 0)
    preambles_delayed1 = [complex(zeros(1,delay1),zeros(1,delay1)), preambles(1,1:end-delay1)];
    preambles_delayed2 = [complex(zeros(1,delay2),zeros(1,delay2)), preambles(1,1:end-delay2)]; 
    y_channel = (h110*preambles + h111*preambles_delayed1 + h112*preambles_delayed2);
else
    y_channel = h110*preambles;
end

preambles_delayed1 = [complex(zeros(1,delay1),zeros(1,delay1)), preambles(1,1:end-delay1)];
preambles_delayed2 = [complex(zeros(1,delay2),zeros(1,delay2)), preambles(1,1:end-delay2)];
y_channel = preambles + preambles_delayed1;%h110*preambles_delayed;

%y_channel = preambles;

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
NIFFT_CORR = 2048;
pdp_freq = ifft(multiplied_sequences,NIFFT_CORR)/Nzc;
pdp_freq_adjusted = pdp_freq;

%% --------------------------------------------------------------------------
if(0) % Apenas um pico por vez
    PDP_ADJUST_FACTOR_0 = (2.292818612477107 - 1.036051241556373i); % delay 0
    PDP_ADJUST_FACTOR_1 = 1; % delay 1
    PDP_ADJUST_FACTOR_2 = 1; % delay 2
    PDP_ADJUST_FACTOR_3 = 1; % delay 3
    PDP_ADJUST_FACTOR_10 = 1; % delay 10
    PDP_ADJUST_FACTOR_13 = 1; % delay 13
    
    pdp_freq_adjusted(710) = pdp_freq_adjusted(710);
    if(delay1 > 0)
        pdp_freq_adjusted(710+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(710+delay_zcorr1);
    end
    
    if(delay2 > 0)
        pdp_freq_adjusted(710+delay_zcorr2) = (PDP_ADJUST_FACTOR_2)*pdp_freq_adjusted(710+delay_zcorr2);
    end
end

%% -------------------------------------------------------------------------
pdp = abs(pdp_freq_adjusted).^2;

stem(0:1:NIFFT_CORR-1,pdp)
%stem(0:1:NIFFT_CORR-1,pdp_freq_adjusted)

h110
h111
h112
pdp_freq_adjusted(1732)
pdp_freq_adjusted(1732+delay_zcorr1)
pdp_freq_adjusted(1732+delay_zcorr2)

