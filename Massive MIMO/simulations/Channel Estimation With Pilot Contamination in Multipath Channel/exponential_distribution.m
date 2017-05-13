clear all;close all;clc

N = 10000000;
x = randn(1,N);
y = randn(1,N);

z = x.^2 + y.^2;

hist(z,1000)