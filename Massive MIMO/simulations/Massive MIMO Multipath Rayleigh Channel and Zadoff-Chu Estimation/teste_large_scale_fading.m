clear all;clc

%% ----------------------------------- Parameters ------------------------------
plot_fig = false;                                                   % Disable/enable figure plotting.

L = 1;                                                              % Number of cells.
K = 5;                                                              % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
NFFT = 2048;                                                        % Number of points used by the OFDM.
modOrd = 2;                                                         % Constellation size = 2^modOrd.
numSym = K*NFFT;                                                    % Number of symbols, i.e., number of terminals.
NCP = 128;                                                          % Number of samples used to create a Extended Cyclic Prefix (CP) according to 3GPP' LTE standards. 12 OFDM symbols per subframe, i.e., 1 ms.
Ts = 1/(15000*NFFT);                                                % System Sampleing Rate.
numSymbInSubframe = 12;                                             % Number of symbols in a subframe. 12 for extended CP.

EbNoVec = -20:1:0;                                                  % Eb/No in dB.
EsN0dB = EbNoVec + 10*log10(NFFT/(NFFT+NCP)) + 10*log10(modOrd);    % converting to symbol to noise ratio
snr = EsN0dB - 10*log10((NFFT/(NFFT+NCP)));                         % Calculate SNR from EsNo in dB.

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

nTotalOfBits = 1e7;
nErrors = 100000;
debug = false;

% Large scale fading.
sshadow = 3;                    % Shadow-fading standard deviation in dB.
gamma = 2.8;                    % Decay exponent: Urban area cellular radio 2.7 - 3.5

% Small scale fading.
delay = [0 0.4]*1e-6;           % Atraso em microsegundos.
gain  = [0 0];                  % Ganho em db
numPaths = length(delay);       % Number of paths per channel.
totalNumPaths = M*K*numPaths;   % Total number of paths between the various antennas and Base Stations.

%% ----------------------------------- Setup Channel -----------------------------------
% Vetor de ganhos.
pos =  round(delay/Ts)+1;               % Posicao dos taps efetivos no vetor.
g   = zeros(1, round(delay(end)/Ts)+1); % +1 serve para incluir o delay 0 no vetor.
for n = 1:length(delay)
    g( pos(n) ) = 10^( gain(n)/10 );
end

fc = 1e9;               % Carrier Freq. in MHz.
c = 3e8;                % Light speed in m/s.
v = 30;                 % Speed in m/s.
Fs = 30.72e6;           % Sampling Freq. in MHz.
Ts = 1/Fs;              % Sampling period in seconds.

fd = (v*fc)/c;          % Doppler frequency.

Pd = 0;                 % Potencia relativa em dB.

% Parameters used for generating the Mulipath Masive MIMO Channel.
Fs_chann = 1000;                    % Channel Sampling Rate in Hz.
Ts_chann = 1/Fs_chann;
N_chann = 256;                      % Number of samples used to sample the channel.
delta_f_chann = Fs_chann/N_chann;   % in Hz.
f = -Fs_chann/2:delta_f_chann:Fs_chann/2;
idx = find(f<fd);
f = f(idx);
idx = find(f>-fd);
f = f(idx);
f = f.';

LS = length(f);
S = 1/pi/fd./sqrt(1 - (f/fd).^2) * 10^(Pd/10);
S = S * LS / sum(S);                % Normalizacao de energia
S1 = S;
S = [S((LS+1)/2:LS); zeros(N_chann-LS,1); S(1:(LS-1)/2)];

% **************** Generate M-MIMO Channel. ****************
rng(55);
x = [(randn((LS-1)/2+1,totalNumPaths,'double') + 1i*randn((LS-1)/2+1,totalNumPaths,'double')); zeros(N_chann-LS,totalNumPaths); (randn((LS-1)/2,totalNumPaths,'double')) + 1i*randn((LS-1)/2,totalNumPaths,'double')]/sqrt(2);
ch = ifft(x .* repmat(sqrt(S),1,totalNumPaths)) * N_chann / sqrt(LS);

% **************** Cell Layout. ****************
% Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
radius = cellHole + (cellRadius-cellHole).*rand(1,K);

% Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
angle = 2*pi*rand(1,K);

% Plot the position of each terminal inside the cell.
figure;
plot(0, 0, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]);
hold on;

Circlex = cellRadius*cos(2*pi*(0:99)/100);
Circley = cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

Circlex = cellHole*cos(2*pi*(0:99)/100);
Circley = cellHole*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'r');

UEx = radius.*cos(angle);
UEy = radius.*sin(angle);
plot(UEx, UEy, 'b*');

grid on;
hold off;

% Large scale fading according to described in Marzetta.
if(1)
    shadowing_att = lognrnd(0,sshadow,1,K);
    path_loss = radius.^gamma;
    largeScaleFading = shadowing_att./path_loss;
end

% Large scale fading according to described in MIMO-OFDM Wireless Communications with MATLAB.
if(0)
    d0 = 100;
    sigma = 2;
    n = 2;
    PL = PL_logdist_or_norm(fc,radius,d0,n,sigma);
    largeScaleFading = 10.^(-0.1*PL);
end

largeScaleFading_sqrt = repmat( sqrt(largeScaleFading), M, 1);

% Apply power delay profile to the channel paths.
for idx_ch = 1 : N_chann
    
    % Atualizacao das matrizes de canal para a celula alvo (l = 1): existe uma matriz de canal por percurso do terminal.
    G = reshape(ch(idx_ch,:), M, K, length(pos));
    G(:,:,1) = (g(pos(1)) * G(:,:,1)) .* largeScaleFading_sqrt;
    for k = 2:length(pos)
        G(:,:,k) = (g(pos(k)) * G(:,:,k)) .* largeScaleFading_sqrt;
    end
    ch(idx_ch,:) = reshape(G, 1, totalNumPaths);
    
    
    DiagMatrix1 = diag(abs((G(:,:,1)' * G(:,:,1))/M));
    DiagMatrix2 = diag(abs((G(:,:,2)' * G(:,:,2))/M));
    diagVector1(idx_ch,:) = (DiagMatrix1);
    diagVector2(idx_ch,:) = (DiagMatrix2);
    
end


%sum(diagVector)

% G = reshape(ch(idx_ch,:), M, K, length(pos));
% DiagMatrix = (abs((G(:,:,1)' * G(:,:,1))/M));
% diagVector = diag(DiagMatrix);

