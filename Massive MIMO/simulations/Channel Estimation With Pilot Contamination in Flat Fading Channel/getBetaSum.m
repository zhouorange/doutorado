function [beta_sum_ik, beta11k] = getBetaSum(beta_i, k_idx, L)

beta_sum_ik = 0;
for l_idx=1:1:L
    beta_sum_ik = beta_sum_ik + beta_i(l_idx,k_idx);
    if(l_idx==1)
        beta11k = beta_i(1,k_idx);
    end
end
