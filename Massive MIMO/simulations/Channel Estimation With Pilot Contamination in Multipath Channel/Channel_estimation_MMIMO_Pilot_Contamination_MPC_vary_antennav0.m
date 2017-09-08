clear all;close all;clc

rng(1)

SNR = 0;          % Signal-to-noise ratio in dB.

M = 10:10:140;   % Number of antennas.
K = 10;            % Number of single-antenna users.
P = 20;            % Channel Length (Finite Impulse Response - FIR).
L = 7;             % Number of cells.

N = getPilotLength(K,P);  % Pilot length is set according to K and P.
q = 10.^(SNR./10);            % Uplink pilot power.

NUM_ITER = 10000;

isFixedValue = true;
if(isFixedValue)
    a = 0.05;
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
set(gca,'yscale','log','xlim',[M(1) M(end)],'ylim',[0.22 0.36]);
xlabel('M'); ylabel('MSE'); set(fdee_figure,'NumberTitle','off');
set(fdee_figure, 'renderer', 'zbuffer');
titleStr = sprintf('SNR: %d - K: %d - P: %d - L: %d - N: %d',SNR,K,P,L,N);
title(titleStr);

ls_error_vec = zeros(1,length(M));
mmse_error_vec = zeros(1,length(M));
prop_1_error_vec = zeros(1,length(M));
prop_2_error_vec = zeros(1,length(M));
prop_3_error_vec = zeros(1,length(M));
prop_4_error_vec = zeros(1,length(M));
for m_idx=1:1:length(M)
    
    epsilon_hat = 0;
    for q_idx=1:1:length(q)
        
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
        
        ls_error = complex(zeros(P,P),zeros(P,P));
        mmse_error = complex(zeros(P,P),zeros(P,P));
        prop_error1 = complex(zeros(P,P),zeros(P,P));
        prop_error2 = complex(zeros(P,P),zeros(P,P));
        prop_error3 = complex(zeros(P,P),zeros(P,P));
        prop_error4 = complex(zeros(P,P),zeros(P,P));
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
            
            % 3) Proposed Estimation 1.
            % Z11_prop = (M*beta111*(Z11_ls * ((aux)^(-1)))); % Ideia original: resultado ruim.
            ep1 = aux(P,P); % Primeira ideia: resultado bom, bate com teórico.
            Z11_prop1 = (M(m_idx)*beta111*Z11_ls)./ep1;
            
            prop_error1 = prop_error1 + ((Z11_prop1'-C111')*(Z11_prop1-C111));
            
            epsilon_hat = epsilon_hat + (ep1/M(m_idx));
            
            % 4) Proposed estimation 2.
            ep2 = (sum(diag(aux))./P); % Segunda ideia: resultado melhor que o teórico.
            Z11_prop2 = (M(m_idx)*beta111*Z11_ls)./ep2;
            
            prop_error2 = prop_error2 + ((Z11_prop2'-C111')*(Z11_prop2-C111));
            
            % 5) Proposed estimation 3.
            ep3 = diag(diag(aux));
            inv_ep3 = (ep3)^(-1);
            Z11_prop3 = (M(m_idx)*beta111*Z11_ls)*inv_ep3;
            
            prop_error3 = prop_error3 + ((Z11_prop3'-C111')*(Z11_prop3-C111));
            
            % 6) Proposed estimation 4.
            num_diags = 4;
            aux_diag = diag(aux);
            ep4 = (sum(aux_diag(1:num_diags))./num_diags);
            Z11_prop4 = (M(m_idx)*beta111*Z11_ls)./ep4;
            
            prop_error4 = prop_error4 + ((Z11_prop4'-C111')*(Z11_prop4-C111));
        end
        
        epsilon_hat = epsilon_hat/NUM_ITER;
        
        fprintf('SNR: %d - M: %d\n',SNR(q_idx), M(m_idx));
        
        % Least Squares (LS) Simulation Error. (Monte Carlo)
        ls_error = ls_error./(M(m_idx) * NUM_ITER);
        ls_error_vec(m_idx) = sum(diag(ls_error))/P;
        
        % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
        mmse_error = mmse_error./(M(m_idx) * NUM_ITER);
        mmse_error_vec(m_idx) = sum(diag(mmse_error))/P;
        
        % Proposed Estimator 1 Simulation Error. (Monte Carlo)
        prop_error1 = prop_error1./(M(m_idx) * NUM_ITER);
        prop_1_error_vec(m_idx) = sum(diag(prop_error1))/P;
        
        % Proposed Estimator 2 Simulation Error. (Monte Carlo)
        prop_error2 = prop_error2./(M(m_idx) * NUM_ITER);
        prop_2_error_vec(m_idx) = sum(diag(prop_error2))/P;
        
        % Proposed Estimator 3 Simulation Error. (Monte Carlo)
        prop_error3 = prop_error3./(M(m_idx) * NUM_ITER);
        prop_3_error_vec(m_idx) = sum(diag(prop_error3))/P;
        
        % Proposed Estimator 4 Simulation Error. (Monte Carlo)
        prop_error4 = prop_error4./(M(m_idx) * NUM_ITER);
        prop_4_error_vec(m_idx) = sum(diag(prop_error4))/P;
        
    end
    
    % Plot results.
    semilogy(M(1:m_idx),real(mmse_error_vec(1:m_idx)),'r*', ...
        M(1:m_idx),real(ls_error_vec(1:m_idx)),'b*', ...
        M(1:m_idx),real(prop_1_error_vec(1:m_idx)),'k*', ...
        M(1:m_idx),real(prop_2_error_vec(1:m_idx)),'kv', ...
        M(1:m_idx),real(prop_3_error_vec(1:m_idx)),'ko', ...
        M(1:m_idx),real(prop_4_error_vec(1:m_idx)),'ks' ,'MarkerSize',7)
    legend('MMSE (sim)', 'LS (sim)', 'Prop. 1 (sim)', 'Prop. 2 (sim)', 'Prop. 3 (sim)', 'Prop. 4 (sim)');
    drawnow;
    
end

semilogy(M,real(mmse_error_vec),'r-', ...
    M,real(ls_error_vec),'b-', ...
    M,real(prop_1_error_vec),'k-', ...
    M,real(prop_2_error_vec),'k-', ...
    M,real(prop_3_error_vec),'k-', ...
    M,real(prop_4_error_vec),'k-' ,'MarkerSize',7)
hold off;

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v0_%s.fig',K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v0_%s.mat',K,P,L,N,timeStamp);
save(fileName);
