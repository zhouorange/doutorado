clear all;close all;clc

rng(1)

SNR = 10;           % Signal-to-noise ratio in dB.

M = 90;             % Number of antennas.
K = 10;             % Number of single-antenna users.
L = 7;              % Number of cells.

N = K;              % Pilot length is set according to K and P.
q = 10.^(SNR./10);  % Uplink pilot power.

a = logspace(-3, 1, 40);        % Interferance leval value.
beta111 = 1;

NUM_ITER = 10000;

% Generate pilot signals.
S = generatePilotMatrixFFT(N,K);

theoretical_ls_error = zeros(1,length(a));
theoretical_mmse_error = zeros(1,length(a));
theoretical_proposed_error = zeros(1,length(a));
theoretical_proposed_approx_error = zeros(1,length(a));
ls_error_vec = zeros(1,length(a));
mmse_error_vec = zeros(1,length(a));
prop_error_vec = zeros(1,length(a));
for a_idx=1:1:length(a)    
    ls_error = 0;
    mmse_error = 0;
    prop_error = 0;
    for numIter = 1:1:NUM_ITER
        
        beta_sum = 0;
        sum_G = zeros(M,N);
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
        W1 = (1/sqrt(2))*complex(randn(M,N),randn(M,N));
        
        % received pilot symbols at BS.
        Y1 = sum_G + W1;
        
        % ******* LS ********
        Z11_ls = (1/(sqrt(q)*N))*Y1*S(:,1);       
        ls_error = ls_error + ((Z11_ls'-Gil(:,1,1)')*(Z11_ls-Gil(:,1,1)));
        
        % ******* MMSE ********
        Z11_mmse = (beta111/epsilon11)*Z11_ls;
        mmse_error = mmse_error + ((Z11_mmse'-Gil(:,1,1)')*(Z11_mmse-Gil(:,1,1)));
        
        % ******* Prop. ********
        Z11_prop = M*beta111*Z11_ls/((Z11_ls')*(Z11_ls));
        prop_error = prop_error + ((Z11_prop'-Gil(:,1,1)')*(Z11_prop-Gil(:,1,1)));
    end
    
    fprintf('a: %d\n',a(a_idx));
    
    % Least Squares (LS) Simulation Error. (Monte Carlo)
    ls_error = ls_error./(M * NUM_ITER);
    ls_error_vec(a_idx) = ls_error;
    
    % Least Squares (LS) Theoretical error.
    theoretical_ls_error(a_idx) = epsilon11 - beta111;
    
    % Minimum Mean Squared Error (MMSE) Simulation Error. (Monte Carlo)
    mmse_error = mmse_error./(M * NUM_ITER);
    mmse_error_vec(a_idx) = mmse_error;
    
    % Minimum Mean Squared Error (MMSE) Theoretical error.
    theoretical_mmse_error(a_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
    
    % Proposed estimator Simulation Error (Monte Carlo)
    prop_error = prop_error./(M * NUM_ITER);
    prop_error_vec(a_idx) = prop_error;   
    
    % Proposed estimator Theoretical error.
    theta_ik = calculateTheta_ikv2(beta111, epsilon11, M);
    theoretical_proposed_error(a_idx) = ((M/(M-1))*((beta111^2)/epsilon11)) + beta111 - 2*beta111*theta_ik;
    
    % Proposed estimator Approximated error.
    theoretical_proposed_approx_error(a_idx) = (beta111*((((2-M)*beta111)/((M-1)*epsilon11)) + 1));
end

fontSize = 10;

fdee_figure = figure;
loglog(a,real(theoretical_mmse_error),'-b');
hold on;
loglog(a,real(mmse_error_vec),'b*','MarkerSize',7);
loglog(a,real(theoretical_ls_error),'-g','MarkerSize',7);
loglog(a,real(ls_error_vec),'g*','MarkerSize',7);
loglog(a,real(theoretical_proposed_error),'-r','MarkerSize',7);
loglog(a,real(theoretical_proposed_approx_error),'rx','MarkerSize',7);
loglog(a,real(prop_error_vec),'ro','MarkerSize',7);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('MMSE (ana)','MMSE (sim)','LS (ana)', 'LS (sim)', 'Prop. (ana)', 'Prop. (approx.)', 'Prop. (sim)', 'Location','northwest');
strText = sprintf('M = %d, q = %1.0f dB',M,q);
x1 = a(1)+(30/100)*a(1);
y1 = 1.5;
text(x1,y1,strText,'FontSize', fontSize)

% Get timestamp for saving files.
timeStamp = datestr(now,30);
fileName = sprintf('Figure3_MSE_versus_Cross_Cell_Interference_Level_M%d_K%d_L%d_N%d_v0_%s.fig',M,K,L,N,timeStamp);
savefig(fdee_figure,fileName);

% Save workspace to .mat file.
clear fdee_figure
fileName = sprintf('Figure3_MSE_versus_Cross_Cell_Interference_Level_M%d_K%d_L%d_N%d_v0_%s.mat',M,K,L,N,timeStamp);
save(fileName);
