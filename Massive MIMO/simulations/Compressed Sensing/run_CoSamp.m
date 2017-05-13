close all;clear all;clc

idx = 0;
error_cse_cosamp = Inf;
while(error_cse_cosamp>1e-100)
    
    rng('shuffle');
    
    [error_cse_cosamp, min_error, ppos_min] = teste_ompv6_pilot_dist_cosampv1();
    fprintf(1,'****************************************** min error: %d\n',min_error);
    
    idx = idx + 1;
end


