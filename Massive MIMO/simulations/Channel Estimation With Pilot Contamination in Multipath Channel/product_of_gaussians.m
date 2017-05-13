clear all;close all;clc

N = 10; %how many pairs of variables
maxMean = 5; %i'll make N variables with different means; this is the max
meansX = maxMean*rand(1,N); %make E[X_n]
variancesX = rand(1,N); %Var[X_n]
meansY = maxMean*rand(1,N); %make E[Y_n]
variancesY = rand(1,N); %Var[Y_n]
Ntrials = 1000000; %now many trials of: X_1*Y_1+X_2*Y_2+...+X_N*Y_N
X = randn(Ntrials,N);Y = randn(Ntrials,N);%make normal(0,1); scaled below
for i = 1:N %run through the N pairs and set the means and variances
    X(:,i) = sqrt(variancesX(i))*X(:,i)+meansX(i);
    Y(:,i) = sqrt(variancesY(i))*Y(:,i)+meansY(i);
end
Z = sum(X.*Y,2);hist(Z,100);
title(sprintf('Mean of %d trials = %5.3f; Sum of means of each distr = %5.3f',Ntrials,sum(meansX.*meansY),mean(Z)));
sum(meansX.^2.*variancesY+meansY.^2.*variancesX+variancesY.*variancesY)
var(Z)