function r = channel_calc(tau_max,s,h)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%     Time-Variant Channel with additive Gaussian noise (channel realizations are calculated)
%%%
%%%%%%%%%%%%%%%%%%%
%%%
%%%     r           ...     received signal (column vector) 
%%%     h           ...     channel realizations ( = time dependent filters)
%%%                       size(h) = [Td, tau_max] (Td denotes duration of
%%%                       transmission and is calculated from tau_max and length of s)
%%%     tau_max     ...     maximum delay
%%%     s           ...     transmit signal (column vector) 
%%%
%%%
%%%     we implement: 
%%%
%%%     r(m)= sum_{l=0}^{tau_max-1} h(m,l) s(m-l)
%%%
%%%     which reads in MATLAB, where vector indices start at 1:
%%%     r(m)= sum_{l=1}^{tau_max} h(m,l) s(m-l+1)

Td=length(s);

s_zeros = [zeros(tau_max-1,1) ; s ; zeros(tau_max-1,1)];
r = zeros(Td, 1);
for m = 1:Td 
    indx = m + tau_max - (1:tau_max);
    r(m) = h(m, :) * s_zeros(indx);
end

