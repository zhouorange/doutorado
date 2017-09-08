function [S] = generatePilotMatrixFFT(K, k_idx)

S = fft(eye(K));

S = S(:,k_idx);

