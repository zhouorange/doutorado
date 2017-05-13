clear all;clc;close all;

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
v = [10 11 12];
Ncs = 13;
K = 3;
NFRAME = 30720;

signal = complex(0,0)*zeros(1,NIFFT);

preambles = complex(zeros(K,NFRAME),zeros(K,NFRAME));

show_figures = false;

tx_shift_adjust = 10;

%% ****************************** PRACH Transmission ******************************
for i=1:1:K
    % ------- Generate local Zadoff-Chu sequence -------
    n = [0:1:(Nzc-1)];
    xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);
    
    Cv = v(i)*Ncs;
    xuv = xu_root(mod((n+Cv),Nzc)+1);    
    
    % ------- Generate base-band signal -------
    Xuv = fft(xuv,Nzc);
    signal(1:Nzc) = Xuv;
    bb_signal = ifft(signal, NIFFT);
    
    % ------- Add CP -------
    % Preamble format 0:
    Nseq = 2048;
    Ncp = 264;
    bb_signal_cp = [bb_signal(NIFFT-Ncp+1:NIFFT), bb_signal, zeros(1,tx_shift_adjust)];
    
    if(show_figures)
        figure(1)
        stem(abs(fft(bb_signal_cp,NIFFT)));
        title('Base-band signal with CP - 2.56 Mbps');
    end
    
    % ------- Up-sampling by a factor of 12 -------
    y = upsample(bb_signal_cp,12);
    
    Hd = butterworth3; % Best response curve: flat!
    y_filtered = filter(Hd,y);
    
    if(show_figures)
        figure(2)
        stem(abs(fft(y,NIDFT)));
        title('Upsampled signal - 2.56 Mbps * 12 = 30.72 Mbps');
        
        figure(3)
        stem(abs(fft(y_filtered((Ncp*12)+1:NIDFT+12*Ncp),NIDFT)));
        title('Low-pass filtered Upsampled signal');
    end
    
    y_filtered_adjusted = y_filtered(12*tx_shift_adjust+1:length(y_filtered));
    
    % ------- Time-domain frequency shift -------
    theta = 7;
    K = 12;
    Nrbsc = 12;
    Nulrb = 25;
    nraprb = 4;
    ko = nraprb*Nrbsc - (Nulrb*Nrbsc)/2;
    Ncp = 3168;
    fo = theta+K*(ko+1/2);
    m = 0:1:NIDFT+Ncp-1;
    k = m-Ncp;
    time_freq_shift = exp((1i*2*pi*fo*k)/NIDFT);
    y_shifted = y_filtered_adjusted.*time_freq_shift;
    
    % ------- Add Guard-band interval -------
    NG = 2976;
    y_shifted = [y_shifted, zeros(1,NG)];
    
    if(show_figures)
        figure(4)
        stem(abs(fft(y_shifted(Ncp+1:NIDFT+Ncp),NIDFT)));
        title('Time Domain Frequency Shifted Base-band Signal');
    end
    
    preambles(i,:) = y_shifted;
end

%% *************************** Rayleigh Channel ***************************

h11 = 0.5 + 1i*0.8;
h12 = 0.2 - 1i*0.3;
h13 = 0.15 + 1i*0.4;

y_channel = h11*preambles(1,:) + h12*preambles(2,:) + h13*preambles(3,:);

%% ****************************** PRACH Reception ******************************

% ------- Remove CP and GB -------
rec_signal = y_channel(Ncp+1:NIDFT+Ncp);

if(show_figures)
    figure(5)
    stem(abs(fft(rec_signal,NIDFT)));
    title('Received base-band signal');
end

% ------- Time-domain frequency shift -------
m = 0:1:NIDFT-1;
time_freq_shift = exp((-1i*2*pi*fo*m)/NIDFT);

rec_signal_shifted = rec_signal.*time_freq_shift;

% Adjust delay added by the downsampling filter.
rx_shift_adjust = 3;
adjusted_rec_signal_shifted = [rec_signal_shifted zeros(1,12*rx_shift_adjust)];

if(show_figures)
    figure(6)
    stem(abs(fft(adjusted_rec_signal_shifted,NIDFT)));
    title('Time-domain frequency shifted Received base-band signal');
end

% ------- Downsampling with polyphase filter -------
downsample_factor = 12;
num = [0.00017939033531018 0.000373331595313973 0.00071274312540008 0.00121789473997835 0.00189381011836236 0.00275388375215944 0.00376067623652285 0.00487234517759502 0.00598957449393566 0.00700728758190838 0.00777151701050356 0.00814045427108632 0.00795185498237852 0.00709158862282701 0.00546778367396934 0.00307861305755767 -1.63753592195977e-005 -0.00363923795041085 -0.00753728454262852 -0.0113456884758382 -0.0146531315425441 -0.0169868383944301 -0.017897105615824 -0.0169574750468651 -0.0138553470030869 -0.00839044785522955 -0.000552643816493893 0.00950073408600361 0.0213970980368926 0.0346051144947382 0.0484238439582872 0.0620729050177951 0.0747101454691587 0.0855405957225286 0.0938411757641215 0.0990638802182588 0.100842037329038 0.0990638802182588 0.0938411757641215 0.0855405957225286 0.0747101454691587 0.0620729050177951 0.0484238439582872 0.0346051144947382 0.0213970980368926 0.00950073408600361 -0.000552643816493893 -0.00839044785522955 -0.0138553470030869 -0.0169574750468651 -0.017897105615824 -0.0169868383944301 -0.0146531315425441 -0.0113456884758382 -0.00753728454262852 -0.00363923795041085 -1.63753592195977e-005 0.00307861305755767 0.00546778367396934 0.00709158862282701 0.00795185498237852 0.00814045427108632 0.00777151701050356 0.00700728758190838 0.00598957449393566 0.00487234517759502 0.00376067623652285 0.00275388375215944 0.00189381011836236 0.00121789473997835 0.00071274312540008 0.000373331595313973 0.00017939033531018];
hm = mfilt.firdecim(downsample_factor,num);

downsampled_signal = filter(hm,adjusted_rec_signal_shifted);

adjusted_downsampled_signal = downsampled_signal(rx_shift_adjust+1:length(downsampled_signal));

if(show_figures)
    figure(8)
    stem(abs(fft(adjusted_downsampled_signal,NIFFT)));
    title('Downsampled received signal');
end

% ------- FFT received signal -------
rec_fft = fft(adjusted_downsampled_signal,NIFFT); % The FFT result is scaled by NIFFT. It is also done in the system generator model.

% ------- Sub-carrier de-mapping -------
rec_Xuv = rec_fft(1:Nzc);

% ------- Generate local Zadoff-Chu root sequence with fixing cyclic shift due to the delay caused by the filter(s) -------
Xu_root = fft(xu_root, Nzc);

% ------- Multiply Local Zadoff-Chu root sequence by received sequence -------
conj_Xu_root = conj(Xu_root);
multiplied_sequences = (rec_Xuv.*conj_Xu_root);

% ------- Squared modulus used for peak detection -------
PDP_ADJUST_FACTOR = (10.606925276014076 - 1i*6.875162692444854);
pdp_freq = ifft(multiplied_sequences,Nzc)/Nzc;
pdp_freq_adjusted = PDP_ADJUST_FACTOR*pdp_freq;
pdp = abs(pdp_freq_adjusted).^2;

stem(0:1:Nzc-1,pdp)
