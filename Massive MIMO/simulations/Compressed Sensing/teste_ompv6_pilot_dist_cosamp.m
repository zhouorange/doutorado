clear all;close all;clc;

N = 2048;               % number of OFDM subcarriers.
NCP = 128;              % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
Fs = 30.72e6;
T = 1/Fs;

%% ********* Pilot sequence. *********
% Comment: The vector P addresses pilot sequences.
M = 9; %9                                         % Number of pilots. Length of the training sequence.
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

ppos_min = zeros(1,M);
ppos = 0;
min_error = Inf;
error_cse_cosamp = Inf;
while(error_cse_cosamp>1e-100)
    
    %% ******** Pilot Positions *********
    while(length(unique(ppos)) ~= M)
        ppos = randi(N,1,M);
    end
    ppos = sort(ppos);
    
    %ppos = [506 797 833 855 932 993 1144 1366 1550 1607 1661 1770 1809 1872 1913 1991 2024]; % best positions for M=17
    %ppos = [41 75 153 241 411 542 594 745 1298 1324 1375 1440 1588 1787 1800 1926]; % best positions for M=16
    %ppos = [32 386 420 445 465 513 889 934 1021 1071 1094 1204 1610 1744 1791]; % best positions for M=15
    %ppos = [246 391 492 659 661 848 1406 1462 1491 1564 1604 1774 1924 1951]; % best positions for M=14
    %ppos = [324 426 436 550 929 948 1189 1299 1314 1368 1577 1742 1841]; %best positions for M=13
    %ppos = [77 84 128 565 1076 1167 1297 1324 1413 1425 1488 1777]; %best positions for M=12
    %ppos = [305 334 495 518 580 766 1137 1445 1779 1830 1943]; %best positions for M=11
    %ppos = [396 503 957 1049 1134 1338 1477 1516 1585 1941 2013]; %best positions for M=11
    %ppos = [396 503 957 1049 1338 1477 1516 1585 1941 2013]; %best positions for M=10
    ppos = [74 96 499 717 877 1566 1674 1710 1866]; %best positions for M=9
    
    %% ********* Fourier Basis *********
    Fl = F(ppos,:);
    
    % Comments: The channel "h" has a structured model. You can insert more realistics channel models.
    H = Fl*h; % H is the FFT of the channel H computed at the pilot's frequencies.
    
    %% ************* Transmission ****************
    
    % Transmitted signal in frequency domain.
    s = Pdiag*H;
    
    % received signal computed at pilot's frequencies.
    y = s;
    
    % Compressed Sensing
    A = Pdiag*Fl;
    [g_hat_cse_cosamp] = CoSaMP(A,y,K);
    
    % Error
    error_cse_cosamp = norm(h(1:length(g_hat_cse_cosamp))-g_hat_cse_cosamp)^2/norm(h(1:length(g_hat_cse_cosamp)))^2;
    
    if(error_cse_cosamp < min_error)
        min_error = error_cse_cosamp;
        ppos_min = ppos;
    end
    
    fprintf(1,'error: %d - min error: %d\n',error_cse_cosamp,min_error);
a=1;
end