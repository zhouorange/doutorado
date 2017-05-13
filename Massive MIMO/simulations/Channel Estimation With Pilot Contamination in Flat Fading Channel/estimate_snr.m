clear all;close all;clc

rng(1)

SNR = 10;           % Signal-to-noise ratio in dB.

M = 30;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.
q = 10.^(SNR./10);  % Uplink pilot power.

a = 0.05; %logspace(-3, 1, 40);        % Interferance leval value.
beta111 = 1;

NUM_ITER = 100000;

% Generate pilot signals.
S = sqrt(N)*eye(N);%generatePilotMatrixFFT(N);

for iter=1:1:10
    snr_est = 0;
    sum_G_mean = 0;
    W1_mean = 0;
    for a_idx=1:1:length(a)
        for numIter = 1:1:NUM_ITER
            
            beta_sum = 0;
            sum_G = zeros(M,K);
            Gil = zeros(M,K,L);
            for l=1:1:L
                
                % Generate channels.
                beta = a(a_idx);
                if(l == 1)
                    beta = 1;
                end
                betaMatrix = sqrt(beta)*eye(K);
                Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
                
                % Summation of all channels.
                sum_G = sum_G + Gil(:,:,l)*(S');
                
                % Summation of betas.
                beta_sum = beta_sum + beta;
                
            end
            
            % Factor.
            epsilon11 = (beta_sum + 1/(q*N));
            
            % Apply squared pilot power.
            sum_G = sqrt(q)*sum_G;
            
            % Generate noise.
            W1 = (1/sqrt(2))*complex(randn(M,K),randn(M,K));
            
            % received pilot symbols at BS.
            Y1 = sum_G + W1;
            
            sum_G_mean = sum_G_mean + (sum_G*sum_G'); %(sum_G'*sum_G);
            W1_mean = W1_mean + (W1*W1'); %(W1'*W1);
        end
    end
    sum_G_mean = sum_G_mean/NUM_ITER;
    W1_mean = W1_mean/NUM_ITER;
    snr_est = sum_G_mean/W1_mean;
    fprintf('a: %d - Estimated SNR: %1.4f\n',a(a_idx),trace(snr_est)/M);
end
