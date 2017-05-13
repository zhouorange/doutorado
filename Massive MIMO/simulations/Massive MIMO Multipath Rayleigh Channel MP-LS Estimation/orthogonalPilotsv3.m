clear all;close all;clc

Np = 128;

% BPSK pilots.
ip = rand(1,Np) > 0.5;
P = 2*ip-1;

%P = [1 -1 -1 -1 -1 -1 1 1 1 -1 -1 -1 1 -1 -1 -1 1 1 1 -1 1 -1 1 -1 1 -1 -1 1 1 -1 1 1];

n = 0:1:Np-1;
xuv = zeros(10,Np);
for Cv=0:1:Np-1
    xuv(Cv+1,:) = P(mod((n-Cv),Np)+1);
end

res = (xuv*xuv')/Np;

res2 = abs(sum((xuv(1,:).*conj(xuv(2,:)))))/Np;

%res1 = abs(sum((F(10,:).*conj(F(11,:)))))/NFFT;

a=1;