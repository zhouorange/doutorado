function [S] = generatePilotMatrixFFTv3(K)

S = (1/sqrt(K))*fft(eye(K));

