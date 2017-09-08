function sum_of_all_betas = getSumOfAllBetas(beta,L,K)

sum_of_all_betas = 0;
for l_idx=1:1:L
    for k_idx=1:1:K
        sum_of_all_betas = sum_of_all_betas + beta(l_idx,k_idx);
    end
end