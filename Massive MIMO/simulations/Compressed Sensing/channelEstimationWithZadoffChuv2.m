clear all;clc;close all;

% ------- Definitions -------
u = 25;
Nzc = 63;
NIFFT = 2048;
NIDFT = 24576;
NSUBFRAME = 30720;
Ncp = 3168;
NG = 2976;
v = [0 1 2];
Ncs = 13;
position = mod(Nzc-v.*Ncs,Nzc) + 1; % position of the start of the window. Plus one to correct address matlab vectors.

K = 1; % Number of Tranmitters, i.e., UEs or number of single-antenna terminals.

prach_offset = 10;

signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);

preambles = complex(zeros(K,NSUBFRAME),zeros(K,NSUBFRAME));

show_figures = false;

% ------- Generate Root Zadoff-Chu sequence. -------
xu_root = complex(zeros(1,Nzc),zeros(1,Nzc));

n = [0:1:30];
xu_root(1,1:31) = exp(-1i*(pi*u.*n.*(n+1))./Nzc);

n = [31:1:62];
xu_root(1,32:63) = exp(-1i*(pi*u.*n.*(n+1).*(n+2))./Nzc);

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
delay_zcorr1 = 600;
delay1 = delay_zcorr1*1;

h110 = 0.5 + 1i*0.8;
h111 = 0.2 + 1i*0.2;

if(delay1 > 0)
    delayed_preambles = h111*[complex(zeros(1,delay1),zeros(1,delay1)), preambles(1,1:end-delay1)];
else
    delayed_preambles = complex(zeros(1,NSUBFRAME),zeros(1,NSUBFRAME));
end
direct_preambles = h110*preambles(1,:);
y_channel = direct_preambles + delayed_preambles;

snr = 0;
rx = awgn(y_channel, snr, 0);

rx = y_channel;

%% ****************************** PRACH Reception ******************************

% ------- CP and GT Removal. -------
rec_signal = rx(Ncp+1:NIDFT+Ncp);

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
NIFFT_CORR = Nzc;
zu = ifft(multiplied_sequences,NIFFT_CORR)/Nzc;

stem(0:1:NIFFT_CORR-1,abs(zu).^2)

% pdp_freq_adjusted = zu;

% %% --------------------------------------------------------------------------
% if(0) % Apenas um pico por vez
%     PDP_ADJUST_FACTOR_0 = 1;%0.977180705809313 + 0.003166526587357i; % delay 0    
%     PDP_ADJUST_FACTOR_1 = 0.989338652615384 + 0.152095818347174i; % delay 1   
%     PDP_ADJUST_FACTOR_2 = 0.957498745916891 + 0.301528270811912i; % delay 2
%     PDP_ADJUST_FACTOR_3 = 0.904911629047376 + 0.445660126688108i; % delay 3
%     PDP_ADJUST_FACTOR_10 = 0.050047412688511 + 1.101839958815881i; % delay 10
%     PDP_ADJUST_FACTOR_13 = -0.473961715008805 + 1.083869828182508i; % delay 13
%     
%     pdp_freq_adjusted(position(1)) = PDP_ADJUST_FACTOR_0*pdp_freq_adjusted(position(1));
%     pdp_freq_adjusted(position(2)) = PDP_ADJUST_FACTOR_0*pdp_freq_adjusted(position(2));
%     pdp_freq_adjusted(position(3)) = PDP_ADJUST_FACTOR_0*pdp_freq_adjusted(position(3));
%     if(delay1 > 0)
%         pdp_freq_adjusted(position(1)+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(position(1)+delay_zcorr1);
%         pdp_freq_adjusted(position(2)+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(position(2)+delay_zcorr1);
%         pdp_freq_adjusted(position(3)+delay_zcorr1) = (PDP_ADJUST_FACTOR_1)*pdp_freq_adjusted(position(3)+delay_zcorr1);
%     end
% end
% 
% %% -------------------------------------------------------------------------
% pdp = abs(pdp_freq_adjusted).^2;
% 
% figure;
% stem(0:1:NIFFT_CORR-1,pdp)
% %stem(0:1:NIFFT_CORR-1,pdp_freq_adjusted)
% 
% figure;
% pdp_freq_adjusted_inv=1./pdp_freq_adjusted;
% stem(0:1:NIFFT_CORR-1,abs(pdp_freq_adjusted_inv))
% 
% % fprintf(1,'h110: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h110),position(1),num2str(pdp_freq_adjusted(position(1))),abs(h110-pdp_freq_adjusted(position(1))));
% % fprintf(1,'h120: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h120),position(2),num2str(pdp_freq_adjusted(position(2))),abs(h120-pdp_freq_adjusted(position(2))));
% % fprintf(1,'h130: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h130),position(3),num2str(pdp_freq_adjusted(position(3))),abs(h130-pdp_freq_adjusted(position(3))));
% % 
% % fprintf(1,'h111: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h111),(position(1)+delay_zcorr1),num2str(pdp_freq_adjusted(position(1)+delay_zcorr1)),abs(h111-pdp_freq_adjusted(position(1)+delay_zcorr1)));
% % fprintf(1,'h121: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h121),(position(2)+delay_zcorr1),num2str(pdp_freq_adjusted(position(2)+delay_zcorr1)),abs(h121-pdp_freq_adjusted(position(2)+delay_zcorr1)));
% % fprintf(1,'h131: %s - pdp_freq_adjusted(%d): %s - error: %f\n',num2str(h131),(position(3)+delay_zcorr1),num2str(pdp_freq_adjusted(position(3)+delay_zcorr1)),abs(h131-pdp_freq_adjusted(position(3)+delay_zcorr1)));
