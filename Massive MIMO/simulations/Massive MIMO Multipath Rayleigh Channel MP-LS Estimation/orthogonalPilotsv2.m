clear all;close all;clc

NFFT = 32;

F = fft(eye(NFFT));

P = F(3,:);

n = 0:1:NFFT-1;
xuv = zeros(10,NFFT);
for Cv=0:1:NFFT-1
    xuv(Cv+1,:) = P(mod((n-Cv),NFFT)+1);
end

%res = xuv*xuv';

res2 = abs(sum((xuv(1,:).*conj(xuv(15,:)))))/NFFT;

%res1 = abs(sum((F(10,:).*conj(F(11,:)))))/NFFT;

a=1;