clear all;close all;clc

rng(1)

SNR = -10:2:22;          % Signal-to-noise ratio in dB.

M = 30;%30; %85;        % Number of antennas.
K = 10;             % Number of single-antenna users.
P = 20;              % Channel Length (Finite Impulse Response - FIR).
L = 7;              % Number of cells.

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

ls_error_vec = zeros(1,length(q));
theoretical_ls_error = zeros(1,length(q));
mmse_error_vec = zeros(1,length(q));
theoretical_mmse_error = zeros(1,length(q));
prop_1_error_vec = zeros(1,length(q));
prop_2_error_vec = zeros(1,length(q));
prop_3_error_vec = zeros(1,length(q));
prop_4_error_vec = zeros(1,length(q));
theoretical_proposed_error = zeros(1,length(q));
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
        
        % 3) Proposed Estimation 1.
        % Z11_prop = (M*beta111*(Z11_ls * ((aux)^(-1)))); % Ideia original: resultado ruim.
        ep1 = aux(P,P); % Primeira ideia: resultado bom, bate com teórico.
        Z11_prop1 = (M*beta111*Z11_ls)./ep1;
        
        prop_error1 = prop_error1 + ((Z11_prop1'-C111')*(Z11_prop1-C111));
        
        epsilon_hat = epsilon_hat + (ep1/M);
        
        % 4) Proposed estimation 2.
        ep2 = (sum(diag(aux))./P); % Segunda ideia: resultado melhor que o teórico.
        Z11_prop2 = (M*beta111*Z11_ls)./ep2;
        
        prop_error2 = prop_error2 + ((Z11_prop2'-C111')*(Z11_prop2-C111));
        
        % 5) Proposed estimation 3.
        ep3 = diag(diag(aux));
        inv_ep3 = (ep3)^(-1);
        Z11_prop3 = (M*beta111*Z11_ls)*inv_ep3;
        
        prop_error3 = prop_error3 + ((Z11_prop3'-C111')*(Z11_prop3-C111));
            
        % 6) Proposed estimation 4.
        num_diags = 4;
        aux_diag = diag(aux);
        ep4 = (sum(aux_diag(1:num_diags))./num_diags);
        Z11_prop4 = (M*beta111*Z11_ls)./ep4;
        
        prop_error4 = prop_error4 + ((Z11_prop4'-C111')*(Z11_prop4-C111));
    end
    
    epsilon_hat = epsilon_hat/NUM_ITER;
    
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
    prop_error1 = prop_error1./(M * NUM_ITER);
    prop_1_error_vec(q_idx) = sum(diag(prop_error1))/P;
    
    theta_ik = calculateTheta_ikv1(beta111, epsilon11, M);
    %theta_ik = 0.354626695759462; % para valores de M > 85
    theoretical_proposed_error(q_idx) = ((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
    
    % Proposed Estimator 2 Simulation Error. (Monte Carlo)
    prop_error2 = prop_error2./(M * NUM_ITER);
    prop_2_error_vec(q_idx) = sum(diag(prop_error2))/P;
    
    % Proposed Estimator 3 Simulation Error. (Monte Carlo)
    prop_error3 = prop_error3./(M * NUM_ITER);
    prop_3_error_vec(q_idx) = sum(diag(prop_error3))/P;
    
    % Proposed Estimator 4 Simulation Error. (Monte Carlo)
    prop_error4 = prop_error4./(M * NUM_ITER);
    prop_4_error_vec(q_idx) = sum(diag(prop_error4))/P;
    
end

fdee_figure = figure;
semilogy(SNR,theoretical_mmse_error,'--r');
hold on;
semilogy(SNR,real(mmse_error_vec),'r*','MarkerSize',7);
semilogy(SNR,theoretical_ls_error,'--b','MarkerSize',7);
semilogy(SNR,real(ls_error_vec),'b*','MarkerSize',7);
semilogy(SNR,theoretical_proposed_error,'--k','MarkerSize',7);
semilogy(SNR,real(prop_1_error_vec),'k*','MarkerSize',7);
semilogy(SNR,real(prop_2_error_vec),'-kv','MarkerSize',7);
semilogy(SNR,real(prop_3_error_vec),'-ko','MarkerSize',7);
semilogy(SNR,real(prop_4_error_vec),'-ks','MarkerSize',7);
hold off
grid on;
%axis([-10 24 0.23 0.55])
titleStr = sprintf('M: %d - K: %d - P: %d - L: %d - N: %d',M,K,P,L,N);
title(titleStr);
axis([-10 22 0.2 1.304])
xlabel('SNR [dB]')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)','Prop. (ana)','Prop. 1 (sim)','Prop. 2 (sim)','Prop. 3 (sim)', 'Prop. 4 (sim)');

if(0)
    figure;
    semilogy(SNR,theoretical_mmse_error,'r');
    hold on;
    semilogy(SNR,theoretical_ls_error,'-bv','MarkerSize',7);
    semilogy(SNR,theoretical_proposed_error,'-ko','MarkerSize',7);
    semilogy(SNR,prop_error_vec,'kx','MarkerSize',7);
    hold off
    grid on;
    axis([-10 24 0.2 1.07])
    xlabel('SNR [dB]')
    ylabel('MSE')
    legend('MMSE','LS','Prop. (ana)','Prop. (sim)');
end

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v9_%s.fig',M,K,P,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('channel_estimation_mse_vs_tx_snr_M%d_K%d_P%d_L%d_N%d_v9_%s.mat',M,K,P,L,N,timeStamp);
save(fileName);
