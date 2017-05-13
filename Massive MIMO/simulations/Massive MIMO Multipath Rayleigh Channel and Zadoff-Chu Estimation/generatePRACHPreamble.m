function [preambles] = generatePRACHPreamble(NT)

% ----------------------------- Definitions -------------------------------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
Ncs = 13;
signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);
y_shifted = complex(0,0)*zeros(1,30720);
preambles = complex(0,0)*zeros(NT,30720);

show_figures = false;

tx_shift_adjust = 10;

% ------------------- Generate local Zadoff-Chu sequence ------------------
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);

for nt = 0:1:(NT-1)
        
    Cv = nt*Ncs;
    xuv = xu_root(mod((n+Cv),Nzc)+1);

    % ----------------------- Generate base-band signal -----------------------
    Xuv = fft(xuv,Nzc);
    signal(1 : Nzc) = Xuv;
    bb_signal = ifft(signal, NIFFT);
    
    % ---------------------------------- Add CP -------------------------------
    % Preamble format 0:
    Ncp = 264;
    bb_signal_cp = [bb_signal(NIFFT-Ncp+1:NIFFT), bb_signal, zeros(1,tx_shift_adjust)];
    
    % --------------------- Up-sampling by a factor of 12 ---------------------
    y = upsample(bb_signal_cp,12);
    
    Hd = butterworth3; % Best response curve: flat!
    y_filtered = filter(Hd,y);
    y_filtered_adjusted = y_filtered(12*tx_shift_adjust+1:length(y_filtered));
    
    % ---------------------- Time-domain frequency shift ----------------------
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
    
    % ------------------------- Add Guard-band interval -----------------------
    NG = 2976;
    y_shifted = [y_shifted, zeros(1,NG)];
    
    preambles(nt+1,:) = y_shifted;

end