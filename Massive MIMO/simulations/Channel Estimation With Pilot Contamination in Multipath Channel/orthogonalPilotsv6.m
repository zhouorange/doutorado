clear all;close all;clc

Np = 32;

F = fft(eye(Np));

P1 = F(3,:);
P2 = F(7,:); 
P3 = F(17,:); 

n = 0:1:Np-1;
xuv1 = zeros(Np,Np);
xuv2 = zeros(Np,Np);
xuv3 = zeros(Np,Np);
for Cv=0:1:Np-1
    xuv1(Cv+1,:) = P1(mod((n-Cv),Np)+1);
    xuv2(Cv+1,:) = P2(mod((n-Cv),Np)+1);
    xuv3(Cv+1,:) = P3(mod((n-Cv),Np)+1);
end

XUV = [xuv1 xuv2 xuv3];

res = (XUV')*xuv2;


res2 = sum(P2.*conj(P1));

res3 = abs(res)/Np;