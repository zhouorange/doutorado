clear all;close all;clc

% Construct the compressed sensing problem
n = 200;
m = 100;
A = randn(m,n);
u = full(sprandn(n,1,.1));
b = A*u + 0.01*randn(m,1); % add noise

epsilon = 0.02;

x=sdpvar(n,1);
solvesdp(norm(A*x-b)<=epsilon,norm(x,1));