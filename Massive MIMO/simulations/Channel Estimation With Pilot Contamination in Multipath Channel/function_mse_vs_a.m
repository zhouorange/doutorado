function[theoretical_mmse_error, theoretical_ls_error, prop_4_ana_error_vec] = function_mse_vs_a(SNR, M, K, P, L, a)

rng(1)

N = getPilotLength(K,P);            % Pilot length is set according to K and P.
q = 10.^(SNR./10);                  % Uplink pilot power.

isFixedValue = true;

theoretical_ls_error = zeros(1,length(a));
theoretical_mmse_error = zeros(1,length(a));
prop_4_ana_error_vec = zeros(1,length(a));
for a_idx=1:1:length(a)
    
    % Calculate Beta.
    [beta_sum, beta111, beta] = generateBetav0(a(a_idx), L, K, isFixedValue);
    
    for q_idx=1:1:length(q)
        
        epsilon11 = (beta_sum + 1/(q(q_idx)*N));
                
        fprintf('a: %d\n',a(a_idx));
        
        % Least Squares (LS) Theoretical error.
        theoretical_ls_error(a_idx) = epsilon11 - beta111;
        
        % Minimum Mean Squared Error (MMSE) Theoretical error.
        theoretical_mmse_error(a_idx) = (beta111/epsilon11) * (epsilon11 - beta111);
        
        % Proposed Estimator Theoretical error.
        prop_4_ana_error_vec(a_idx) = beta111*( ( ((2-M*P)*beta111)/((P*M-1)*epsilon11)  ) + 1);
        
    end
end
