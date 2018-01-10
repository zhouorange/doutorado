clear all;close all;clc

M = 10;
a = 2;

% numIter = 10000000;
% z1_acc = zeros(1,numIter);
% z2_acc = zeros(1,numIter);
% zM_acc = zeros(1,numIter);
% for iter=1:1:numIter
%     
%     x = sqrt(a)*(1/sqrt(2))*complex(randn(M,1),randn(M,1));
%     
%     denominator = (x'*x);
%     z1 = (x(1,1) / denominator);
%     z2 = (x(2,1) / denominator);
%     zM = (x(M,1) / denominator);
%     
%     z1_acc(iter) = z1;
%     z2_acc(iter) = z2;
%     zM_acc(iter) = zM;
%     
% end
% 
% z1_real_mean = mean(real(z1_acc));
% z1_imag_mean = mean(imag(z1_acc));
% z1_var = var(z1_acc);
% 
% z2_real_mean = mean(real(z2_acc));
% z2_imag_mean = mean(imag(z2_acc));
% z2_var = var(z2_acc);
% 
% zM_real_mean = mean(real(zM_acc));
% zM_imag_mean = mean(imag(zM_acc));
% zM_var = var(zM_acc);

n = M;
u = -0.4:0.01:0.4;
fu = zeros(1,length(u));
for u_idx=1:1:length(u)
    bessel1 = besseli(n,(1./(4*(u(u_idx).^2))));
    bessel2 = besseli(n+1,(1./(4*(u(u_idx).^2))));
    fu(u_idx) = exp(-1./(4*(u(u_idx).^2)))*( ( (8*n*u(u_idx).^2 - 1)*bessel1 + bessel2 ) / (4*(abs(u(u_idx)).^3)) );
end
u = 0.1;
fu_ = (1./(2*sqrt(2*pi)))*(4*n-2+(-8*(n.^3)+12*(n.^2)-3)*(u.^2)+(16*(n.^5)-40*(n.^3)+9*(n))*(u.^4));

figure(1);
%histogram(real(z1_acc),'Normalization','pdf');
hold on
plot(u,fu,'LineWidth',2)
hold off
