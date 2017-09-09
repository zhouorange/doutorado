clear all;close all;clc

sigma = linspace(0.81e-9,200e-9,10000); %rms delay spread

Bc = 1./(5.*(sigma));


semilogy(sigma,Bc,'LineWidth',1)
hold on
semilogy(sigma,(1.08e6)*ones(1,length(sigma)),'LineWidth',1)
semilogy((1/(5*1.08e6))*ones(1,length(sigma)),linspace(Bc(1),Bc(length(Bc)),10000),'LineWidth',1)
xlabel('RMS delay spread [s]');
ylabel('Bc [Hz]');
axis([sigma(1), sigma(length(sigma)), 1e6 2.47e8])
grid on
hold off