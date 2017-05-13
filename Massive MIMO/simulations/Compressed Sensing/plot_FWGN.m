% plot_FWGN.m
clear, clf
fm=100; scale=1e-6; % Maximum Doppler frequency and mu
ts_mu=50; ts=ts_mu*scale; fs=1/ts; % Sampling time/frequency
Nd=1e6; % Number of samples
% obtain the complex fading channel
[h,Nfft,Nifft,doppler_coeff] = FWGN_model(fm,fs,Nd);
subplot(211), plot([1:Nd]*ts,10*log10(abs(h)))
str=sprintf('Clarke/Gan Model, f_m=%d[Hz], T_s=%d[us]',fm,ts_mu);
title(str), axis([0 0.5 -30 5])
subplot(223), hist(abs(h),50), subplot(224), hist(angle(h),50)
mean(h)