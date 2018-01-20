clear all;close all;clc

M = 100;
zeta = 2;

beta_iik = 1;

numSamples = 1000000;

zik = zeros(M,numSamples);
giik = zeros(M,numSamples);
gcovmatrix = zeros(M,M);
for i=1:1:numSamples
    zik(:,i) = sqrt(zeta)*(1/sqrt(2))*complex(randn(1,M),randn(1,M));
    
    norm2 = (zik(:,i)')*zik(:,i);
    giik(:,i) = beta_iik*M*( zik(:,i) / norm2 );
    
    gcovmatrix = gcovmatrix + (giik(:,i)*giik(:,i)');
    
end

gcovmatrix = gcovmatrix/numSamples;