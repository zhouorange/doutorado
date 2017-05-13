clear all;close all;clc

rng(1)

SNR = 20;          % Signal-to-noise ratio in dB.

M = 30;            % Number of antennas.
K = 10;            % Number of single-antenna users.
P = 1:2:56;        % Channel Length (Finite Impulse Response - FIR).
L = 7;             % Number of cells.

q = 10.^(SNR./10); % Uplink pilot power.

a = 0.05;

NUM_ITER = 10000;

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
xlabel('P'); ylabel('MSE'); set(fdee_figure,'NumberTitle','off');
set(fdee_figure, 'renderer', 'zbuffer');
titleStr = sprintf('SNR: %d - K: %d - M: %d - L: %d',SNR,K,M,L);
title(titleStr);

ls_error_vec = zeros(1,length(P));
mmse_error_vec = zeros(1,length(P));
prop_error_vec = zeros(1,length(P));
N_vector = zeros(1,length(P));
for p_idx=1:1:length(P)
    
    N = getPilotLength(K,P(p_idx));  % Pilot length is set according to K and P.
    N_vector(p_idx) = N;
    
    for m_idx=1:1:length(M)
        
        for q_idx=1:1:length(q)
            
            epsilon11 = (beta_sum + 1/(q(q_idx)*N));
            
            ls_error = complex(zeros(P(p_idx),P(p_idx)),zeros(P(p_idx),P(p_idx)));
            mmse_error = complex(zeros(P(p_idx),P(p_idx)),zeros(P(p_idx),P(p_idx)));
            prop_error = complex(zeros(P(p_idx),P(p_idx)),zeros(P(p_idx),P(p_idx)));
            for iter=1:1:NUM_ITER
                
                % Pilot signal generation
                S = [];
                for k_idx=1:1:K
                    s_aux = generatePilotMatrixv1(N, P(p_idx), k_idx);
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
                        g = (1/sqrt(2))*complex(randn(M(m_idx),P(p_idx)),randn(M(m_idx),P(p_idx)));
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
                
                % 3) Proposed Estimation: P = P(p_idx).
                num_diags = P(p_idx);
                aux_diag = diag(aux);
                ep1 = (sum(aux_diag(1:num_diags))./num_diags);
                Z11_prop1 = (M(m_idx)*beta111*Z11_ls)./ep1;
                
                prop_error = prop_error + ((Z11_prop1'-C111')*(Z11_prop1-C111));
            end
            
            fprintf('SNR: %d - P: %d\n',SNR(q_idx), P(p_idx));
            
            % Minimum Mean Squared Error (MMSE) Theoretical error.
            mmse_error_vec(p_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
            
            % Proposed Estimator 1 Simulation Error. (Monte Carlo)
            prop_error = prop_error./(M(m_idx) * NUM_ITER);
            prop_error_vec(p_idx) = sum(diag(prop_error))/P(p_idx);
            
        end
        
        % Plot results. 
        semilogy(P(1:p_idx),real(mmse_error_vec(1:p_idx)),'r*', ...
            P(1:p_idx),real(prop_error_vec(1:p_idx)),'k*' ,'MarkerSize',7)
        legend('MMSE (ideal)', 'Prop. (sim)');
        drawnow;
        
    end
    
end

semilogy(P,real(mmse_error_vec),'r-', ...
    P,real(prop_error_vec),'k-' ,'MarkerSize',7)
hold off;

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_P_K%d_M%d_L%d_N%d_v2_%s.fig',K,M,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_P_K%d_M%d_L%d_N%d_v2_%s.mat',K,M,L,N,timeStamp);
save(fileName);
