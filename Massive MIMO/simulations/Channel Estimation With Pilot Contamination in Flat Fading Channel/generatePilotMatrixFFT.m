function [S] = generatePilotMatrixFFT(N,K)

if(N >= K) 
    eye_size = N;
else
    eye_size = K;
end

S = fft(eye(eye_size));

S = S(1:N,1:K);
