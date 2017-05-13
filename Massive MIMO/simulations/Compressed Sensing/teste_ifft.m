clear all;clc

NFFT = 2048;

IF = ifft(eye(NFFT));

X = rand(1,NFFT);

x1 = ifft(X,NFFT);

x2 = IF*X.';

x2 = x2.';

sum(abs(x1-x2).^2)/NFFT