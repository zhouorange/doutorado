clear all;close all;clc

rng(1)

SNR = 10;                   % Signal-to-noise ratio (SNR) in dB.
linear_SNR = 10.^(SNR./10); % Linear SNR value.

M = 70;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.

a = 0.05;           % Constant beta value.
beta111 = 1;

q = linear_SNR/(N*(1 + a*(L-1)));  % Uplink pilot power.

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

% Simulation loop starts here.
numTrainVectors = 100000;
numTestVectors = 10000;
train_data = zeros(numTrainVectors,M*N*2);
test_data = zeros(numTestVectors,M*N*2);
train_label = zeros(numTrainVectors,2*M);
test_label = zeros(numTestVectors,2*M);
for q_idx=1:1:length(q)
    
    %% Train vectors.
    for trainIter = 1:1:numTrainVectors
        
        beta_sum = 0;
        sum_G = zeros(M,N);
        Gil = zeros(M,K,L);
        % Iterate over all cells (L) in the assumed system.
        % Here we consider the target cell, i, is 1, i.e., i = 1.
        for l=1:1:L
            
            % Generate channel matrix G_{il}.
            beta = a;
            if(l == 1)
                beta = 1;
            end
            betaMatrix = sqrt(beta)*eye(K);
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
            
            % Summation of all channels.
            sum_G = sum_G + Gil(:,:,l)*(S');
            
            % Summation of betas.
            beta_sum = beta_sum + beta;
            
            % This is the label.
            if(l == 1)
                g_111 = Gil(:,1,l);
            end
            
        end
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot-sequence symbols at BS #1, which is the target cell, i.e., i = 1.
        Y1 = sum_G + W1;
        
        idx = 0;
        for y_line_idx=1:1:size(Y1,1)
            for y_col_idx=1:1:size(Y1,2)
                idx = idx + 1;
                train_data(trainIter,idx) = real(Y1(y_line_idx,y_col_idx));
                idx = idx + 1;
                train_data(trainIter,idx) = imag(Y1(y_line_idx,y_col_idx));
            end
        end
        
        idx = 0;
        for g_line_idx=1:1:size(g_111,1)
            for g_col_idx=1:1:size(g_111,2)
                idx = idx + 1;
                train_label(trainIter,idx) = real(g_111(g_line_idx,g_col_idx));
                idx = idx + 1;
                train_label(trainIter,idx) = imag(g_111(g_line_idx,g_col_idx));        
            end
        end
        
    end
    
    %% Test vectors.
    for testIter = 1:1:numTestVectors
        
        beta_sum = 0;
        sum_G = zeros(M,N);
        Gil = zeros(M,K,L);
        % Iterate over all cells (L) in the assumed system.
        % Here we consider the target cell, i, is 1, i.e., i = 1.
        for l=1:1:L
            
            % Generate channel matrix G_{il}.
            beta = a;
            if(l == 1)
                beta = 1;
            end
            betaMatrix = sqrt(beta)*eye(K);
            Gil(:,:,l) = (1/sqrt(2))*complex(randn(M,K),randn(M,K))*betaMatrix;
            
            % Summation of all channels.
            sum_G = sum_G + Gil(:,:,l)*(S');
            
            % Summation of betas.
            beta_sum = beta_sum + beta;
            
            % This is the label.
            if(l == 1)
                g_111 = Gil(:,1,l);
            end
            
        end
        
        % Factor.
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot-sequence symbols at BS #1, which is the target cell, i.e., i = 1.
        Y1 = sum_G + W1;
        
        idx = 0;
        for y_line_idx=1:1:size(Y1,1)
            for y_col_idx=1:1:size(Y1,2)
                idx = idx + 1;
                test_data(testIter,idx) = real(Y1(y_line_idx,y_col_idx));
                idx = idx + 1;
                test_data(testIter,idx) = imag(Y1(y_line_idx,y_col_idx));
            end
        end
        
        idx = 0;
        for g_line_idx=1:1:size(g_111,1)
            for g_col_idx=1:1:size(g_111,2)
                idx = idx + 1;
                test_label(testIter,idx) = real(g_111(g_line_idx,g_col_idx));
                idx = idx + 1;
                test_label(testIter,idx) = imag(g_111(g_line_idx,g_col_idx));        
            end
        end
        
    end
    
end

fileName = sprintf('data_set_M_%d_K_%d_SNR_%d_static_scenario_1.mat',M,K,SNR);
save(fileName,'train_data','train_label','test_data','test_label','-v7')

