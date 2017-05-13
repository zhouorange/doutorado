clear all;clc;

size_s = 10000000;

% noise
rng('shuffle');
sigma = 0.1;



%noise = sqrt(sigma/2)*(randn(size(s))+1i*randn(size(s)));

SNR = 3; % SNR given in dB.
linearSNR = 10^(-SNR/20);
%noise = (randn(1,size_s) + 1i*randn(1,size_s)) / sqrt(2) / ( linearSNR );

%noise = (randn(1,size_s) + 1i*randn(1,size_s)) / sqrt(2);

noise = linearSNR*((randn(1,size_s) + 1i*randn(1,size_s)) / sqrt(2));

sum(abs(noise).^2)/size_s