clear all;close all;clc;

N = 2048;               % number of OFDM subcarriers.
NCP = 128;              % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
Fs = 30.72e6;
T = 1/Fs;

%% ********* Pilot sequence. *********
% Comment: The vector P addresses pilot sequences.
M = 32;                                         % Number of pilots. Length of the training sequence.
%ppos = ceil(linspace(1,N,M));                  % Position of the pilot in the frequency domain.
%ppos = 1:N/M:N;

fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:M-1)/M);               % Normalized pilot signal, i.e., unit power.
%P = 1/sqrt(2)*(randn(1,M) + 1i*randn(1,M));

Pdiag = diag(P);

F = fft(eye(N));

%% ********* Channel ************
h0 = 0.5 + 1i*0.8;
h1 = 0.7 - 1i*0.3;
h2 = 0.1 - 1i*0.4;

energy = (abs(h0).^2 + abs(h1).^2 + abs(h2).^2) / 3;
normalization_factor = 1/sqrt(energy);

h0 = h0*normalization_factor;
h1 = h1*normalization_factor;
h2 = h2*normalization_factor;

channel_gains = [h0, h1, h2];

h = complex(zeros(N,1),zeros(N,1));
pos = [1 24 31];
L = pos(length(pos));
for idx=1:1:length(pos)
    h(pos(idx)) = channel_gains(idx);
end

channel_energy = sum(abs(h).^2)/length(pos);

K = length(pos);                                % Number of non-zero entries

SNR = 20;

ppos_min = zeros(1,M);
min_error = Inf;
error_ce = Inf;
while(error_ce>8e-15)
    
    %% ******** Pilot Positions *********
    ppos = 0;
    while(length(unique(ppos)) ~= M)
        ppos = randi(N,1,M);
    end
    ppos = sort(ppos);
    
    %ppos = [258 305 640 912 1349 1407 1417 1752]; % best pilot positions for M=8. (error: 7.968304e-15)
    %ppos = [525 592 864 1140 1257 1402];% best pilot positions for M=6. (error: 1.158463e-15)
    
    %% ********* Fourier Basis *********
    Fl = F(ppos,:);
    
    % Comments: The channel "h" has a structured model. You can insert more realistics channel models.
    H = Fl*h; % H is the FFT of the channel H computed at the pilot's frequencies.
    
    %% ************* Transmission ****************
    
    % Transmitted signal in frequency domain.
    s = Pdiag*H;
    
    linearSNR = 10^(-SNR/20);
    noise = linearSNR*((randn(size(s)) + 1i*randn(size(s))) / sqrt(2));
    
    % received signal computed at pilot's frequencies.
    y = s + noise;
    
    %% ************ Convex Estimation (CE) ***************
    % Sparse Channel estimation
    epsilon = 0.2;
    A = Pdiag*Fl;
    
    x = sdpvar(N,1);
    solvesdp(norm(A*x-y)<=epsilon,norm(x,1));
    g_hat_ce = double(x);
    
    % Error
    error_ce = norm(h(1:length(g_hat_ce))-g_hat_ce)^2/norm(h(1:length(g_hat_ce)))^2;
    
    if(error_ce < min_error)
        min_error = error_ce;
        ppos_min = ppos;
    end
    
    fprintf(1,'error: %d - min error: %d\n',error_ce,min_error);
    
end