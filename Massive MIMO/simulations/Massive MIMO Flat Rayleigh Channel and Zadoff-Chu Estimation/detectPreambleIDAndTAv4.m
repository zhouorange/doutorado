% rx - received signal, after AWGN and multipath channel.
% M - Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
% K - Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
% numPaths - Number of paths per channel.
% totalNumPaths  - Total number of paths between the various antennas and Base Stations.
% H - Multipath Fading Channel.
% pos - Effective position of taps within the vector.
function [ID, TA, H_estimated] = detectPreambleIDAndTAv4(rx, M, K, numPaths, pos)

%% ---------- PRACH Definitions ----------
u = 129;
Nzc = 839;
NIDFT = 24576;
Ncp = 3168;
v = [0 5 10 15 20 25 30 35 40 45];
Ncs = 13;
position = mod(Nzc-v.*Ncs,Nzc) + 1; % Position of the start of the window. Plus one to correct address matlab vectors.
prach_offset = 10;

PDP_ADJUST_FACTOR = [1, 0.989338652615384 + 0.152095818347174i, 0.957498745916891 + 0.301528270811912i, 0.904911629047376 + 0.445660126688108i, 0.050047412688511 + 1.101839958815881i, -0.473961715008805 + 1.083869828182508i];

delay_zcorr1 = (pos(numPaths) - 1)/30;

%% ------- Generate Root Zadoff-Chu sequence. -------
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);

%% ****************************** PRACH Reception ******************************

% ------- CP and GT Removal. -------
rec_signal = rx(:,Ncp+1:NIDFT+Ncp);

% ------- Apply DFT to received signal. -------
rec_fft = fft(rec_signal,NIDFT,2);

% ------- Sub-carrier de-mapping. -------
rec_Xuv = rec_fft(:,prach_offset+1:prach_offset+Nzc);

% ------- Apply DFT to Root Zadoff-Chu sequence. -------
Xu_root = fft(xu_root, Nzc);

% ------- Multiply Local Zadoff-Chu root sequence by received sequence. -------
conj_Xu_root = conj(Xu_root);
multiplied_sequences = complex(zeros(M,Nzc),zeros(M,Nzc));
for mm=1:1:M
    multiplied_sequences(mm,:) = (rec_Xuv(mm,:).*conj_Xu_root);
end

% ------- Squared modulus used for peak detection. -------
NIFFT_CORR = 839;
pdp_freq = ifft(multiplied_sequences,NIFFT_CORR,2)/Nzc;
pdp_freq_adjusted = pdp_freq;

%% ------------------- Adjust PDP ----------------
for mm=1:1:M
    for idx=1:1:length(position)
        pdp_freq_adjusted(mm,position(idx)) = PDP_ADJUST_FACTOR(1)*pdp_freq_adjusted(mm,position(idx));
        if(numPaths > 0)
            pdp_freq_adjusted(mm,position(idx)+delay_zcorr1) = (PDP_ADJUST_FACTOR(2))*pdp_freq_adjusted(mm,position(idx)+delay_zcorr1);
        end
    end
end

% ----------- Channel estimation. -----------
H_estimated = zeros(M,K,numPaths);
for mm=1:1:M
    for kk=1:1:K
        for idx=1:1:numPaths
            H_estimated(mm,kk,idx) = pdp_freq_adjusted(mm,(position(kk)+(idx-1)));
        end
    end
end
% -----------------------------------------

% ----------------- Retrieve ID and TA ----------------
ID = zeros(1,K);
TA = zeros(K,numPaths);
for kk=1:1:K
    ID(kk) = v(kk);
    for np_idx=1:1:numPaths
        TA(kk,np_idx) = (np_idx-1)*30;
    end
end
% -----------------------------------------
