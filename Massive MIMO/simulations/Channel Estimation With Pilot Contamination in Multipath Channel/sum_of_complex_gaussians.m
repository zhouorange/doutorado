clear all;close all;clc

N = 1e6;

c1 = (sqrt(1))*(1/(sqrt(2)))*complex(randn(1,N),randn(1,N));
c2 = (sqrt(2))*(1/(sqrt(2)))*complex(randn(1,N),randn(1,N));
c3 = (sqrt(3))*(1/(sqrt(2)))*complex(randn(1,N),randn(1,N));

var(c1+c2+c3)
mean(c1+c2+c3)