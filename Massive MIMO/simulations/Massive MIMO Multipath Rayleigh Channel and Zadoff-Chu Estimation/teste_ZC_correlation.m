clear all;clc;close all;

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
v = [0];
Ncs = 13;
signal = complex(0,0)*zeros(1,NIFFT);
xuv = complex(0,0)*zeros(1,Nzc);

show_figures = false;

tx_shift_adjust = 10;

%% *********************** PRACH Transmission *****************************

% ------------------- Generate local Zadoff-Chu sequence ------------------
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);
for i=1:1:length(v)
    Cv = v(i)*Ncs;
    xuv = xuv + xu_root(mod((n+Cv),Nzc)+1);
end


%% *************************** Rayleigh Channel ***************************
h = (randn() + 1i*randn());

y = h*xuv;

%% ************************* PRACH Reception ******************************

%r = cconv(y,conj(xu_root),Nzc)/Nzc; % FUNCIONA!!!!

Xuv = fft(y,Nzc);       % FUNCIONA!!!!
Xu = fft(xu_root,Nzc);
mult = Xuv .* conj(Xu);
r = ifft(mult,Nzc)/Nzc;

pdp = abs(r).^2;

stem(pdp)