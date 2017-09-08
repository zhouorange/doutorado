clear all;close all;clc

rng(1)

SNR = 0;          % Signal-to-noise ratio in dB.

M = 10:10:140;   % Number of antennas.
K = 10;            % Number of single-antenna users.
P = 20;            % Channel Length (Finite Impulse Response - FIR).
L = 7;             % Number of cells.

N = getPilotLength(K,P);  % Pilot length is set according to K and P.
q = 10.^(SNR./10);            % Uplink pilot power.

a = 0.05;

NUM_ITER = 100000;

isFixedValue = true;
if(isFixedValue)
    beta = a*ones(L,K);
    for ll=1:1:L
        for kk=1:1:K
            if(kk==ll)
                beta(ll,kk) = 1;
            end
        end
    end
else
    beta = abs(randn(L,K));
end
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);

% Set up a figure for visualizing BER results.
fdee_figure = figure; grid on; hold on;
%set(gca,'yscale','log','xlim',[M(1) M(end)],'ylim',[0.22 0.36]);
xlabel('M'); ylabel('MSE'); set(fdee_figure,'NumberTitle','off');
set(fdee_figure, 'renderer', 'zbuffer');
titleStr = sprintf('SNR: %d - K: %d - P: %d - L: %d - N: %d',SNR,K,P,L,N);
title(titleStr);

ls_error_vec = zeros(1,length(M));
mmse_error_vec = zeros(1,length(M));
theoretical_ls_error_vec = zeros(1,length(M));
theoretical_mmse_error_vec = zeros(1,length(M));
for m_idx=1:1:length(M)
    
    ls_error = complex(zeros(P,P),zeros(P,P));
    mmse_error = complex(zeros(P,P),zeros(P,P));
    for q_idx=1:1:length(q)
        
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        for iter=1:1:NUM_ITER
            
            % Pilot signal generation
            S = [];
            for k_idx=1:1:K
                s_aux = generatePilotMatrixv1(N, P, k_idx);
                %s_aux = generatePilotMatrixFFT(N, k_idx);
                if(k_idx==1)
                    S1 = s_aux;
                end
                S = [S, s_aux];
            end
            
            % Noise generation
            W = (1/sqrt(2))*complex(randn(M(m_idx),N),randn(M(m_idx),N));
            
            % Received pilot symbols at BS i.
            Y1 = complex(zeros(M(m_idx),N),zeros(M(m_idx),N));
            for l_idx=1:1:L
                
                C = [];
                for k_idx=1:1:K
                    g = (1/sqrt(2))*complex(randn(M(m_idx),P),randn(M(m_idx),P));
                    c = sqrt(beta(l_idx,k_idx))*g;
                    C = [C, c];
                    
                    if(l_idx==1 && k_idx==1)
                        C111 = c;
                    end
                end
                
                Y1 = Y1 + sqrt(q(q_idx))*C*(S');
                
            end
            
            % Add noise to received signal.
            Y1 = Y1 + W;
            
            % 1) Least Squares (LS) Solution (Estimation).
            Z11_ls = (1/(sqrt(q(q_idx))*N))*Y1*S1;
            
            ls_error = ls_error + ((Z11_ls'-C111')*(Z11_ls-C111));
            
            % 2) Minimum Mean Squared Error (MMSE) Solution (Estimation).
            Z11_mmse = (beta111/epsilon11)*Z11_ls;
            
            mmse_error = mmse_error + ((Z11_mmse'-C111')*(Z11_mmse-C111));
            
            % Auxiliar variable.
            aux = ((Z11_ls')*(Z11_ls));
        end
        
        fprintf('SNR: %d - M: %d\n',SNR(q_idx), M(m_idx));
        
        % Least Squares (LS) Theoretical error.
        theoretical_ls_error_vec(m_idx) = epsilon11 - beta111;
        
        % Minimum Mean Squared Error (MMSE) Theoretical error.
        theoretical_mmse_error_vec(m_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
    end
    
    % LS Estimator Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M(m_idx) * NUM_ITER);
    ls_error_vec(m_idx) = sum(diag(ls_error))/P;
    
    % MMSE Estimator Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M(m_idx) * NUM_ITER);
    mmse_error_vec(m_idx) = sum(diag(mmse_error))/P;
    
    % Plot results.
    semilogy(M(1:m_idx),real(theoretical_mmse_error_vec(1:m_idx)),'r-', ...
        M(1:m_idx),real(mmse_error_vec(1:m_idx)),'r*', ...
        M(1:m_idx),real(theoretical_ls_error_vec(1:m_idx)),'k-', ...
        M(1:m_idx),real(ls_error_vec(1:m_idx)),'k*','MarkerSize',7)
    legend('MMSE (ideal)', 'MMSE (sim)', 'LS (ana)', 'LS (sim)');
    drawnow;
    
end

semilogy(M,real(theoretical_mmse_error_vec),'r-', ...
    M,real(mmse_error_vec),'r*', ...
    M,real(theoretical_ls_error_vec),'k-', ...
    M,real(ls_error_vec),'k*','MarkerSize',7)
hold off;

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_est_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v3_%s_simple.fig',K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_est_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v3_%s_simple.mat',K,P,L,N,timeStamp);
save(fileName);
