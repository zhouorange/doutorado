clear all;clc;

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
sc_offset = 0;
v = [10];
%v = [0, 1, 3, 62, 63];
Ncs = 13;
l = 480;
fc = 0.08333;
signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);
Pfa = 0.0001;
Pfd = 0.001;
trials_counter = 0;
false_alarm_counter = 0;
correct_rejection_counter = 0;
detection_counter = 0;
false_rejection_counter = 0;

send_preamble = true;
show_figures = true;
add_noise = true;
add_offset = true;

if(send_preamble)
    nOfTransmittedPreambles = length(v);
else
    nOfTransmittedPreambles = 0;
end

nTrials = 1;

tx_shift_adjust = 10;

error_counter = 0;
error_detected_preamble_v = [];
error_signature = [];
error_TA = [];
error_max_values = [];
error_detected_preambles = [];
error_number_of_detected_preamble = [];

%% ****************************** PRACH Transmission ******************************

% ------- Generate local Zadoff-Chu sequence -------
delay = 4;
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);
if(send_preamble)
    for i=1:1:length(v)
        Cv = v(i)*Ncs;
        xuv = xuv + xu_root(mod((n+Cv),Nzc)+1);
        
        if(v(i)==3)
            xuv = xuv + xu_root(mod((n+Cv-delay),Nzc)+1);
        end
    end
end

% ------- Generate base-band signal -------
Xuv = fft(xuv,Nzc);
signal(1+sc_offset : Nzc+sc_offset) = Xuv;
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
%y = interp(bb_signal_cp,12,l,fc);
y = upsample(bb_signal_cp,12);

%b = fir1(l,fc);
%freqz(b,1,512);
%y_filtered = filter(b,1,y);

%[b,a] = maxflat(10,2,fc);
%y_filtered = filter(b,a,y);
%fvtool(b,a)

Hd = butterworth3; % Best response curve: flat!
%Hd = maximallyflat2;
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

% normalize = true;
% if(normalize)
%     y_shifted = y_shifted/max(abs(y_shifted));
% end

%% ****************************** AWGN and Delay Channel ******************************


if(add_noise)
%     noise_power = -30; % in dBW
%     corrupted_signal = y_shifted + wgn(1,30720,noise_power,'complex');

    snr = -30; % in dB
    corrupted_signal = awgn(y_shifted,snr,'measured');
else
    corrupted_signal = y_shifted;
end

add_offset = false;
if(add_offset)
    offset = 0.3*max(real(corrupted_signal));
    corrupted_signal = corrupted_signal - (offset+0i);
end

add_ta = true;
if(add_ta)
    TA = 0;
    corrupted_signal = [zeros(1,12*TA),corrupted_signal];
end

%% ****************************** PRACH Reception ******************************

% ------- Remove CP and GB -------
rec_signal = corrupted_signal(Ncp+1:NIDFT+Ncp);

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

% ------- Low-pass Filter -------
% Hdd = butterworth3; % Best response curve: flat!
% rec_signal_filtered = filter(Hdd,rec_signal_shifted);
%
% figure(7)
% stem(abs(fft(rec_signal_filtered,NIDFT)));
% title('Low-pass filtered received signal');

% ------- Downsampling -------
%downsampled_signal = downsample(rec_signal_filtered,12);
%downsampled_signal = downsample(rec_signal_shifted,12);

if(show_figures)
    figure(8)
    stem(abs(fft(adjusted_downsampled_signal,NIFFT)));
    title('Downsampled received signal');
end

% ------- FFT received signal -------
rec_fft = fft(adjusted_downsampled_signal,NIFFT)/NIFFT; % The FFT result is scaled by NIFFT. It is also done in the system generator model.

% ------- Estimate the noise reference from the 25 guard subcarriers. -------
nOfRefSamples = 24;
noisy_samples = [rec_fft(Nzc+1:Nzc+12) rec_fft(NIFFT-11:NIFFT)];
%noisy_samples = rec_fft(Nzc+1:Nzc+nOfRefSamples);
%Zref = ((Nzc/NIFFT)^2)*sum(abs(noisy_samples*NIFFT).^2);
Zref = sum(abs(noisy_samples).^2);

% ------- Calculate the scale factor by making use of the Inverse Fisher CDF. -------
alpha = finv(1-Pfa,2,2*nOfRefSamples)/nOfRefSamples;

% ------- Calculate Threshold -------
adjust_factor = 577000;
treshold = alpha*adjust_factor*Zref;

% ------- Sub-carrier de-mapping -------
rec_Xuv = rec_fft(1:Nzc);

% ------- Generate local Zadoff-Chu root sequence with fixing cyclic shift due to the delay caused by the filter(s) -------
max_Xu_root_value = 28.9655;
fix_factor = 0;
xu_root_modified = xu_root(mod((n-fix_factor),Nzc)+1);
Xu_root = fft(xu_root_modified, Nzc)/max_Xu_root_value; % The Root sequence is scaled by max_Xu_root_value so that its maximum value is 1.

% ------- Multiply Local Zadoff-Chu root sequence by received sequence -------
conj_Xu_root = conj(Xu_root);
multiplied_sequences = (rec_Xuv.*conj_Xu_root);

% ------- Squared modulus used for peak detection -------
pdp_freq = ifft(multiplied_sequences,NIFFT)*(max_Xu_root_value*NIFFT);
pdp = abs(pdp_freq).^2;

%if(show_figures)
    stem(0:1:NIFFT-1,pdp)
%end