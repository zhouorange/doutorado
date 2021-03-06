clear all;clc;close all;

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
NSUBFRAME = 30720;
Ncp = 3168;
NG = 2976;
v = [10 20];
Ncs = 13;
position = [710 580];

K = 2; % Number of Tranmitters, i.e., UEs or number of single-antenna terminals.

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
delay_zcorr1 = 1;
delay1 = delay_zcorr1*30;

h110 = 0.5 + 1i*0.8;
h120 = 0.7 + 1i*0.3;

h111 = 0.2 + 1i*0.2;
h121 = 0.4 + 1i*0.6;

energy0 = (abs(h110).^2 + abs(h120).^2) / 2;
normalization_factor0 = 1/sqrt(energy0);

energy1 = (abs(h111).^2 + abs(h121).^2) / 2;
normalization_factor1 = 1/sqrt(energy1);

h110 = h110*normalization_factor0;
h120 = h120*normalization_factor0;

h111 = h111*normalization_factor1;
h121 = h121*normalization_factor1;

energy0 = (abs(h110).^2 + abs(h120).^2) / 2;
energy1 = (abs(h111).^2 + abs(h121).^2) / 2;

if(delay1 > 0)
    delayed_preambles = (h111*[complex(zeros(1,delay1),zeros(1,delay1)), preambles(1,1:end-delay1)]) + (h121*[complex(zeros(1,delay1),zeros(1,delay1)), preambles(2,1:end-delay1)]);  
    los_signal = (h110*preambles(1,:) + h120*preambles(2,:)); % line of sight.
    %delayed_preambles = (1*[complex(zeros(1,delay1),zeros(1,delay1)), preambles(1,1:end-delay1)]) + (1*[complex(zeros(1,delay1),zeros(1,delay1)), preambles(2,1:end-delay1)]);  
    %los_signal = (1*preambles(1,:) + 1*preambles(2,:)); % line of sight.  
    y_channel = los_signal + delayed_preambles;   
else
    y_channel = h110*preambles(1,:) + h120*preambles(2,:);
end

%y_channel = preambles(1,:) + preambles(2,:);

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
NIFFT_CORR = 839;
pdp_freq = ifft(multiplied_sequences,NIFFT_CORR)/Nzc;
pdp_freq_adjusted = pdp_freq;

%% --------------------------------------------------------------------------
if(1) % Apenas um pico por vez
    PDP_ADJUST_FACTOR_0 = 0.977031037663614 + 0.003275798772592i; % delay 0   
    PDP_ADJUST_FACTOR_1 = 0.989338652615384 + 0.152095818347174i; % delay 1
    PDP_ADJUST_FACTOR_2 = 0.957498745916891 + 0.301528270811912i; % delay 2
    PDP_ADJUST_FACTOR_3 = 0.904911629047376 + 0.445660126688108i; % delay 3
    PDP_ADJUST_FACTOR_10 = 0.050047412688511 + 1.101839958815881i; % delay 10
    PDP_ADJUST_FACTOR_13 = -0.473961715008805 + 1.083869828182508i; % delay 13
    
    pdp_freq_adjusted(position(1)) = PDP_ADJUST_FACTOR_0*pdp_freq_adjusted(position(1));
    pdp_freq_adjusted(position(2)) = PDP_ADJUST_FACTOR_0*pdp_freq_adjusted(position(2));
    if(delay1 > 0)
        pdp_freq_adjusted(position(1)+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(position(1)+delay_zcorr1);
        pdp_freq_adjusted(position(2)+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(position(2)+delay_zcorr1);
    end
    
end

%% -------------------------------------------------------------------------
pdp = abs(pdp_freq_adjusted).^2;

figure;
stem(0:1:NIFFT_CORR-1,pdp)
%stem(0:1:NIFFT_CORR-1,pdp_freq_adjusted)

figure;
pdp_freq_adjusted_inv=1./pdp_freq_adjusted;
stem(0:1:NIFFT_CORR-1,abs(pdp_freq_adjusted_inv))
% 
% h110
% h111
% h112
% pdp_freq_adjusted(710)
% pdp_freq_adjusted(710+delay_zcorr1)
% pdp_freq_adjusted(710+delay_zcorr2)

