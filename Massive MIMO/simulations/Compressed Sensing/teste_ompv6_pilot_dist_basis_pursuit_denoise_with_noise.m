clear all;close all;clc;

N = 2048;               % number of OFDM subcarriers.
NCP = 128;              % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
Fs = 30.72e6;
T = 1/Fs;

%% ********* Pilot sequence. *********
% Comment: The vector P addresses pilot sequences.
M = 16;                                         % Number of pilots. Length of the training sequence.
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

SNR = 10;

opts = spgSetParms('verbosity',0);
ppos_min = zeros(1,M);
min_error = Inf;
error_cse = Inf;
ppos = 0;
while(error_cse>1e-10)
    
    %% ******** Pilot Positions *********
    ppos = 0;
    while(length(unique(ppos)) ~= M)
        ppos = randi(N,1,M);
    end
    ppos = sort(ppos);
    
    % ppos = [482 526 579 669 690 738 932 1007 1018 1074 1075 1401 1543 1828 1878 1895]; % ótima distribuição de frequência (subcarriers) para 16 pilotos.
    % ppos = [42 92 121 185 375 531 892 1058 1224 1384 1466 1490 1544 1753 1868 1916];
    % ppos = [254 382 863 1272 1375 1579]; % ótima distribuição de frequência (subcarriers) para 6 pilotos.
    
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
    %y = s + noise;
    
    y = s;
    
    % Compressed Sensing
    A = Pdiag*Fl;
    
    % Reference: https://www.math.ucdavis.edu/~mpf/spgl1/index.html
    % If there was noise inherent in the signal that we measured, it would be better to instead solve the basis pursuit denoise (BPDN) problem
    % minimize ||x||_1  subject to  ||Ax - b||_2 <= sigma,

    sigma = linearSNR; %      % Desired ||Ax - b||_2
    g_hat_cse = spg_bpdn(A, y, sigma, opts);
    
    % Error
    error_cse = norm(h(1:length(g_hat_cse))-g_hat_cse)^2/norm(h(1:length(g_hat_cse)))^2;
    
    if(error_cse < min_error)
        min_error = error_cse;
        ppos_min = ppos;
        g_hat_opt = g_hat_cse;
    end
    
    fprintf(1,'error: %d - min error: %d\n',error_cse,min_error);
    
end