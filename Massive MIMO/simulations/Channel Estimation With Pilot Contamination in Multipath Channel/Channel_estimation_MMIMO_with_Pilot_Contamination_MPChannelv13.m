clear all;close all;clc

rng(1)

SNR = -12:2:20;%-20:2:22;           % Signal-to-noise ratio in dB.

M = 30;%86;%30; %85;                % Number of antennas.
K = 10;                             % Number of single-antenna users.
P = 20;                             % Channel Length (Finite Impulse Response - FIR).
L = 7;                              % Number of cells.

N = getPilotLength(K,P);            % Pilot length is set according to K and P.
q = 10.^(SNR./10);                  % Uplink pilot power.

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

ls_error_vec = zeros(1,length(q));
theoretical_ls_error = zeros(1,length(q));
mmse_error_vec = zeros(1,length(q));
theoretical_mmse_error = zeros(1,length(q));
prop_1_sim_error_vec = zeros(1,length(q));
prop_2_sim_error_vec = zeros(1,length(q));
prop_3_sim_error_vec = zeros(1,length(q));
prop_4_sim_error_vec = zeros(1,length(q));
prop_1_ana_error_vec = zeros(1,length(q));
prop_2_ana_error_vec = zeros(1,length(q));
prop_3_ana_error_vec = zeros(1,length(q));
prop_4_ana_error_vec = zeros(1,length(q));
for q_idx=1:1:length(q)
    
    epsilon11 = (beta_sum + 1/(q(q_idx)*N));
    
    ls_error = complex(zeros(P,P),zeros(P,P));
    mmse_error = complex(zeros(P,P),zeros(P,P));
    prop_sim_error1 = complex(zeros(P,P),zeros(P,P));
    prop_sim_error2 = complex(zeros(P,P),zeros(P,P));
    prop_sim_error3 = complex(zeros(P,P),zeros(P,P));
    prop_sim_error4 = complex(zeros(P,P),zeros(P,P));
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
        W = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % Received pilot symbols at BS i.
        Y1 = complex(zeros(M,N),zeros(M,N));
        for l_idx=1:1:L
            
            C = [];
            for k_idx=1:1:K
                g = (1/sqrt(2))*complex(randn(M,P),randn(M,P));
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
        
        % 3) Proposed Estimation 1: avg P = 1.
        num_diags = 1;
        aux_diag = diag(aux);
        ep1 = (sum(aux_diag(1:num_diags))./num_diags);
        Z11_prop1 = (M*beta111*Z11_ls)./ep1;
        
        prop_sim_error1 = prop_sim_error1 + ((Z11_prop1'-C111')*(Z11_prop1-C111));
        
        % 4) Proposed estimation 2: avg P = 5.
        num_diags = 5;
        aux_diag = diag(aux);
        ep2 = (sum(aux_diag(1:num_diags))./num_diags);
        Z11_prop2 = (M*beta111*Z11_ls)./ep2;
        
        prop_sim_error2 = prop_sim_error2 + ((Z11_prop2'-C111')*(Z11_prop2-C111));
        
        % 5) Proposed estimation 3: avg P = 10.
        num_diags = 10;
        aux_diag = diag(aux);
        ep3 = (sum(aux_diag(1:num_diags))./num_diags);
        Z11_prop3 = (M*beta111*Z11_ls)./ep3;
        
        prop_sim_error3 = prop_sim_error3 + ((Z11_prop3'-C111')*(Z11_prop3-C111));
        
        % 6) Proposed estimation 4: avg P = 20.
        num_diags = 20;
        aux_diag = diag(aux);
        ep4 = (sum(aux_diag(1:num_diags))./num_diags);
        Z11_prop4 = (M*beta111*Z11_ls)./ep4;
        
        prop_sim_error4 = prop_sim_error4 + ((Z11_prop4'-C111')*(Z11_prop4-C111));
    end
    
    fprintf('SNR: %d\n',SNR(q_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M * NUM_ITER);
    ls_error_vec(q_idx) = sum(diag(ls_error))/P;
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error(q_idx) = epsilon11 - beta111;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M * NUM_ITER);
    mmse_error_vec(q_idx) = sum(diag(mmse_error))/P;
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error(q_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
    
    % Proposed Estimator 1 Simulation Error. (Monte Carlo)
    prop_sim_error1 = prop_sim_error1./(M * NUM_ITER);
    prop_1_sim_error_vec(q_idx) = sum(diag(prop_sim_error1))/P;
    
    % Proposed Estimator Theoretical error - P = 1.
    Pavg = 1;
    prop_1_ana_error_vec(q_idx) = beta111*( ( ((2-M*Pavg)*beta111)/((Pavg*M-1)*epsilon11)  ) + 1);
    
    % Proposed Estimator 2 Simulation Error. (Monte Carlo)
    prop_sim_error2 = prop_sim_error2./(M * NUM_ITER);
    prop_2_sim_error_vec(q_idx) = sum(diag(prop_sim_error2))/P;
    
    % Proposed Estimator Theoretical error - P = 5.
    Pavg = 5;
    prop_2_ana_error_vec(q_idx) = beta111*( ( ((2-M*Pavg)*beta111)/((Pavg*M-1)*epsilon11)  ) + 1);    
    
    % Proposed Estimator 3 Simulation Error. (Monte Carlo)
    prop_sim_error3 = prop_sim_error3./(M * NUM_ITER);
    prop_3_sim_error_vec(q_idx) = sum(diag(prop_sim_error3))/P;

    % Proposed Estimator Theoretical error - P = 10.
    Pavg = 10;
    prop_3_ana_error_vec(q_idx) = beta111*( ( ((2-M*Pavg)*beta111)/((Pavg*M-1)*epsilon11)  ) + 1);
    
    % Proposed Estimator 4 Simulation Error. (Monte Carlo)
    prop_sim_error4 = prop_sim_error4./(M * NUM_ITER);
    prop_4_sim_error_vec(q_idx) = sum(diag(prop_sim_error4))/P;
    
    % Proposed Estimator Theoretical error - P = 20.
    Pavg = 20;
    prop_4_ana_error_vec(q_idx) = beta111*( ( ((2-M*Pavg)*beta111)/((Pavg*M-1)*epsilon11)  ) + 1);
    
end

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'--r');
hold on;
semilogy(SNR,real(mmse_error_vec),'r*','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'--b','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b*','MarkerSize',7);
semilogy(SNR,real(prop_1_ana_error_vec),'--k','MarkerSize',7);
semilogy(SNR,real(prop_1_sim_error_vec),'k*','MarkerSize',7);
semilogy(SNR,real(prop_2_ana_error_vec),'--g','MarkerSize',7);
semilogy(SNR,real(prop_2_sim_error_vec),'g*','MarkerSize',7);
semilogy(SNR,real(prop_3_ana_error_vec),'--m','MarkerSize',7);
semilogy(SNR,real(prop_3_sim_error_vec),'m*','MarkerSize',7);
semilogy(SNR,real(prop_4_ana_error_vec),'--c','MarkerSize',7);
semilogy(SNR,real(prop_4_sim_error_vec),'c*','MarkerSize',7);
hold off
grid on;
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d',M,K,P,L,N);
title(titleStr);
%axis([-10 22 0.2 1.304])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)','Prop. 1 avg=1 (ana)','Prop. 1 avg=1 (sim)','Prop. 2 avg=5 (ana)','Prop. 2 avg=5 (sim)','Prop. 3 avg=10 (ana)','Prop. 3 avg=10 (sim)','Prop. 4 avg=20 (ana)','Prop. 4 avg=20 (sim)');

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v13_%s.fig',M,K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v13_%s.mat',M,K,P,L,N,timeStamp);
save(fileName);
