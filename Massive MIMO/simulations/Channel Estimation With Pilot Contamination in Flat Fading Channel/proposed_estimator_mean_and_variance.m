clear all;close all;clc

M = 10;
zeta = 2;
beta_iik = 1;

numIter = 10000000;
z1_acc = zeros(1,numIter);
z2_acc = zeros(1,numIter);
zM_acc = zeros(1,numIter);
for iter=1:1:numIter
    
    x = sqrt(zeta)*(1/sqrt(2))*complex(randn(M,1),randn(M,1));
    
    denominator = (x'*x);
    z1 = M*beta_iik*(x(1,:) / denominator);
    z2 = M*beta_iik*(x(2,:) / denominator);
    zM = M*beta_iik*(x(M,:) / denominator);
    
    z1_acc(iter) = z1;
    z2_acc(iter) = z2;
    zM_acc(iter) = zM;
    
end

z1_real_mean = mean(real(z1_acc));
z1_imag_mean = mean(imag(z1_acc));
z1_var = var(z1_acc);

z2_real_mean = mean(real(z2_acc));
z2_imag_mean = mean(imag(z2_acc));
z2_var = var(z2_acc);

zM_real_mean = mean(real(zM_acc));
zM_imag_mean = mean(imag(zM_acc));
zM_var = var(zM_acc);

theoretical_var = (M/(M-1))*((beta_iik^2)/zeta);

figure(1)
hist(real(z1_acc),1000)
figure(2)
hist(imag(z1_acc),1000)