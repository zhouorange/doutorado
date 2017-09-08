function [S] = generatePilotMatrixFFTv2(K)

S = fft(eye(K));

S = S(2:K,2:K);
