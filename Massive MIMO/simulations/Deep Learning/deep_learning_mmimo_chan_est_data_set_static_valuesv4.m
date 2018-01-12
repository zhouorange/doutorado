clear all;close all;clc

rng(1)

numTrainVectors = 10000;
numTestVectors = 1000;
numPredictionVectors = 100;

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
train_data = zeros(numTrainVectors*N,M*2);
train_label = zeros(numTrainVectors*N,M*2);
test_data = zeros(numTestVectors*N,M*2);
test_label = zeros(numTestVectors*N,M*2);
prediction_data = zeros(numPredictionVectors*N,M*2);
prediction_label = zeros(numPredictionVectors*N,M*2);

for q_idx=1:1:length(q)
    
    %% Train vectors.
    % loop starts here.
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
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot-sequence symbols at BS #1, which is the target cell, i.e., i = 1.
        Y1 = sum_G + W1;
        
        % Generate training data.        
        line_idx = N*(trainIter-1);
        for y_col_idx=1:1:size(Y1,2)
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(Y1,1)
                column_idx = column_idx + 1;
                train_data(line_idx,column_idx) = real(Y1(y_line_idx,y_col_idx));
                column_idx = column_idx + 1;
                train_data(line_idx,column_idx) = imag(Y1(y_line_idx,y_col_idx));
            end
        end
        
        % Generate label data.        
        line_idx = N*(trainIter-1);
        for y_col_idx=1:1:N
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(g_111,1)
                column_idx = column_idx + 1;
                train_label(line_idx,column_idx) = real(g_111(y_line_idx,1));
                column_idx = column_idx + 1;
                train_label(line_idx,column_idx) = imag(g_111(y_line_idx,1));
            end
        end
    end
    
    %% Generate test vectors.
    % loop starts here.
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
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot-sequence symbols at BS #1, which is the target cell, i.e., i = 1.
        Y1 = sum_G + W1;
             
        % Generate training data.
        line_idx = N*(testIter-1);
        for y_col_idx=1:1:size(Y1,2)
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(Y1,1)
                column_idx = column_idx + 1;
                test_data(line_idx,column_idx) = real(Y1(y_line_idx,y_col_idx));
                column_idx = column_idx + 1;
                test_data(line_idx,column_idx) = imag(Y1(y_line_idx,y_col_idx));
            end
        end
        
        % Generate label data.
        line_idx = N*(testIter-1);
        for y_col_idx=1:1:N
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(g_111,1)
                column_idx = column_idx + 1;
                test_label(line_idx,column_idx) = real(g_111(y_line_idx,1));
                column_idx = column_idx + 1;
                test_label(line_idx,column_idx) = imag(g_111(y_line_idx,1));
            end
        end
    end
    
    %% Generate prediction vectors.
    % loop starts here.
    for predictionIter = 1:1:numPredictionVectors
        
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
        
        % Apply squared pilot power.
        sum_G = sqrt(q(q_idx))*sum_G;
        
        % Generate noise.
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot-sequence symbols at BS #1, which is the target cell, i.e., i = 1.
        Y1 = sum_G + W1;
             
        % Generate training data.
        line_idx = N*(predictionIter-1);
        for y_col_idx=1:1:size(Y1,2)
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(Y1,1)
                column_idx = column_idx + 1;
                prediction_data(line_idx,column_idx) = real(Y1(y_line_idx,y_col_idx));
                column_idx = column_idx + 1;
                prediction_data(line_idx,column_idx) = imag(Y1(y_line_idx,y_col_idx));
            end
        end
        
        % Generate label data.
        line_idx = N*(predictionIter-1);
        for y_col_idx=1:1:N
            column_idx = 0;
            line_idx = line_idx + 1;
            for y_line_idx=1:1:size(g_111,1)
                column_idx = column_idx + 1;
                prediction_label(line_idx,column_idx) = real(g_111(y_line_idx,1));
                column_idx = column_idx + 1;
                prediction_label(line_idx,column_idx) = imag(g_111(y_line_idx,1));
            end
        end
    end
    
    %% Save data set for specfic scenario.
    fileName = sprintf('data_set_M_%d_K_%d_SNR_%d_static_scenario_1.mat',M,K,SNR(q_idx));
    save(fileName,'train_data','train_label','test_data','test_label','prediction_data','prediction_label','-v7')
end
