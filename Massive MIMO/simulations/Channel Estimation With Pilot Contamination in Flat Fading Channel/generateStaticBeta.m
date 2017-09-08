function [beta] = generateStaticBeta(a, beta11k, L, K)

beta = zeros(L,K);

for l_idx=1:1:L
    for k_idx=1:1:K   
        if(l_idx==1)
            beta(l_idx,k_idx) = beta11k;
        else
            beta(l_idx,k_idx) = a;
        end        
    end
end
