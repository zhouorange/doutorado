function [beta_sum_ik, beta11k, beta11k_hat1, beta11k_hat2] = getBetaSumv1(L, beta_i, k_idx, var1, var2)

% Beta error.
beta_error_1 = (1 + sqrt(var1)*randn(1));
beta_error_2 = (1 + sqrt(var2)*randn(1));

beta_sum_ik = 0;
for l_idx=1:1:L
    beta_sum_ik = beta_sum_ik + beta_i(l_idx,k_idx);
    if(l_idx==1)
        beta11k = beta_i(1,k_idx);
    end
end

beta11k_hat1 = beta11k*beta_error_1;
beta11k_hat2 = beta11k*beta_error_2;