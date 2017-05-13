clear all;close all;clc;

N = 2048;               % number of OFDM subcarriers.
NCP = 128;              % Number of samples used to create a Cyclic Prefix (CP) according to 3GPP' LTE standards.
Fs = 30.72e6;
T = 1/Fs;

%% ********* Pilot sequence. *********
% Comment: The vector P addresses pilot sequences.
M = 32;                                         % Number of pilots. Length of the training sequence.
% Position of the pilot in the frequency domain.
%ppos = [482 526 579 669 690 738 932 1007 1018 1074 1075 1401 1543 1828 1878 1895]; % ótima distribuição de frequência (subcarriers) para 16 pilotos.
%ppos = [42 92 121 185 375 531 892 1058 1224 1384 1466 1490 1544 1753 1868 1916];
ppos = [186 244 297 311 334 410 417 441 449 760 769 806 813 1001 1027 1030 1048 1105 1146 1176 1203 1350 1380 1470 1604 1606 1610 1624 1731 1898 1954 2016]; % ótima distribuição de frequência (subcarriers) para 32 pilotos.

fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:M-1)/M);               % Normalized pilot signal, i.e., unit power.
%P = 1/sqrt(2)*(randn(1,M) + 1i*randn(1,M));

Pdiag = diag(P);

P = P.';

%% ********* Fourier Basis *********
F = fft(eye(N));        
Fl = F(ppos,:);

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

% Comments: The channel "h" has a structured model. You can insert more realistics channel models.
H = Fl*h; % H is the FFT of the channel H computed at the pilot's frequencies.

% Hl = fft(h,N);
% Hl = Hl(ppos,:);
% H_error = sum(abs(H-Hl).^2)/length(H);

%% ************* Transmission ****************

% Transmitted signal in frequency domain.
s = Pdiag*H;

% noise
%rng(839);

SNR = 2000; % SNR given in dB.
linearSNR = 10^(-SNR/20);
noise = linearSNR*((randn(size(s)) + 1i*randn(size(s))) / sqrt(2));

% received signal computed at pilot's frequencies.
y = s + noise;

%y = s;

% %% ************ Convex Estimation (CE) ***************
% % building dictionary (please check different papers to learn how to build the dictionary)
% GI = NCP;
% TS =  N*T;           % OFDM symbol time (not considering gaurd interval)
% tau_p = linspace(0,GI*T - GI*T./N,N);
% Gamma = exp(-sqrt(-1)*2*pi.*repmat(((1:N).'),1,N)./TS.*repmat(tau_p,N,1));
%     
% % Sparse Channel estimation
% B = 1;
% NP = M;
% Nt = N;
% PP = ppos;
% Pilot = P.';
% pathGains = [h0 h1 h2];
% H_Sparse = zeros(N,B);
% lambda1 = NP*10^(-SNR/10)/sum(abs(pathGains));
% for b = 1 : B
%     A = Gamma(PP,:).*repmat(Pilot(:,b),1,N);
%     cvx_begin quiet
%     variable x(Nt) complex
%     % sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
%     minimize( quad_form(y-A*x,eye(NP))+lambda1*norm(x,1) )
%     cvx_end
%     % building channel at all location (simply from the dictionary)
%     H_Sparse(:,b) = Gamma*x;
% end

%% ************ Compressed Sensing Estimation (CSE) ************

% Estimation of g
%A = Pdiag*Fl(:,1:pos(length(pos)));
%A = Gamma(PP,:).*repmat(Pilot(:,b),1,N);
%A = Pdiag*Fl(:,1:NCP);

alpha = 1;
Nl = alpha*N;
Fll = zeros(length(ppos),NCP*alpha);
for k=1:1:length(ppos)
    Fll(k,:) = exp((-1i*2*pi*(ppos(k)-1)*[0:1:(NCP*alpha-1)])/Nl); 
end

A = Pdiag*Fl;

[g_hat_cse] = OMP_orig(A,y,K);

% Comments: The performance of OMP depends on the choice of F and P. In
% this example, I assume that the channel has Fourier structure. In more
% realistic channel problem, Fourier basis presents a power leakage which
% increases the number of non-zero entries of g. Therefore, the choice of
% the F impacts directly on the algorithm performance.

error_cse = norm(h(1:length(g_hat_cse))-g_hat_cse)^2/norm(h(1:length(g_hat_cse)))^2

%% ************* Least Squares Estimation (LSE) ******************

A = Pdiag*Fl(:,1:pos(length(pos)));
%A = Pdiag*Fl;

invMat = (((A'*A)^(-1))*A');

g_hat_lse = invMat*y;

error_lse = norm(h(1:length(g_hat_lse))-g_hat_lse)^2/norm(h(1:length(g_hat_lse)))^2

%% ************* Mean Squared Error Estimation (MMSE-E) ******************
invMat = (A'*A + ((linearSNR^2)/1)*eye(L))^-1;
invMat = invMat*A';

g_hat_mmsee = invMat*y;

error_mmsee = norm(h(1:length(g_hat_mmsee))-g_hat_mmsee)^2/norm(h(1:length(g_hat_mmsee)))^2


%% ************* Least Squares Estimation (LSEv2) ******************
Hp = y./P;
H = complex(zeros(N,1),zeros(N,1));

H(ppos) = Hp;

H_LSE = ifft(Hp,N);
