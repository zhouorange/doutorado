clear all;close all;clc

NFFT = 2048;                                                        % Number of points used by the OFDM.
Np = 31;
K = 10;
a = (NFFT/sqrt(NFFT-(Np*K)+Np));


Iter=10000;

var_ofdm_rx = 0;
mean_ofdm_rx = 0;
var_r = 0;
mean_r = 0;
for idx=1:1:Iter
    
    % Add channel noise power to faded data.
    r = wgn(1, NFFT, 0);
    
    var_r = var_r + var(r);
    mean_r = mean_r + mean(r);
    
    % Retrieve modulation symbols by applying FFT to received signal.
    b = 1/a;
    ofdm_rx = (a/sqrt(NFFT))*(b*fft(r,NFFT));
    
    %ofdm_rx = (1/sqrt(NFFT))*fft(r,NFFT);
    
    var_ofdm_rx = var_ofdm_rx + var(ofdm_rx);
    mean_ofdm_rx = mean_ofdm_rx + mean(ofdm_rx);
    
end

var_r = var_r/Iter;
mean_r = mean_r/Iter;

var_ofdm_rx = var_ofdm_rx/Iter;
mean_ofdm_rx = mean_ofdm_rx/Iter;

a=1;