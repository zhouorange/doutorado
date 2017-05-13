function [x,l_min,l_max] = MC_demodulator(r,gamma,K,N)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%	Multicarrier Demodulator (Pulse-Shaping)
%%%
%%%%%%%%%%%%%%%%%%%
%%%
%%%   r       ...     received signal (column vector, with length length_r) 
%%%   gamma   ...     receive pulse (vector)
%%%                   assumed to be causal, i.e., 
%%%                   gamma=(gamma[0] ... gamma[length_gamma-1]), otherwise zero
%%%                   in MATLAB - of course - gamma=gamma(1:length_gamma) , i.e., a row
%%%                   vector
%%%   K       ...     number of sub-carriers
%%%   N       ...     discrete symbol period
%%%
%%%   x       ...     matrix consisting of demodulated (data) symbols 
%%%                   number of rows = number of sub-carriers (K)
%%%                   number of columns = number of demodulated (data)
%%%                   symbols = l_max - l_min+1
%%%                   (can be greater than Nsymb, because of fringe effects ("Randeffekte"))   
%%%   l_min   ...     l-index of first column of x (see formula below) (can be negative due to fringe effects)
%%%   l_max   ...     l-index of last column of x (see formula below) (can be greater than Nsymb due to fringe effects) 
%%%
%%% we implement: 
%%%       x_{k,l}=1/sqrt(K) sum_{n=0}^{length_r-1} r(n) gamma*(n-lN) e^{-j 2 pi k/K(n-lN)}


length_gamma=length(gamma);
overlappings=ceil(length_gamma/K);
length_gamma_ext=overlappings*K;
gamma_ext=[gamma zeros(1,(length_gamma_ext-length_gamma))];

r_ext=[zeros((length_gamma_ext-1),1) ; r ; zeros((length_gamma_ext-1),1)];

length_r=length(r);

l_min=ceil(-(length_gamma-1)/N);
l_max=floor((length_r-1)/N);

y=[];
for l=l_min:l_max
    spalt_vec=zeros(K,1);
    for zaehler=0:(overlappings-1)
        spalt_vec=spalt_vec+r_ext((1+length_gamma_ext-1+K*zaehler+l*N):(K+length_gamma_ext-1+K*zaehler+l*N)).*(gamma_ext((1+K*zaehler):(K+K*zaehler)))';
    end
    y=[y spalt_vec];
end   
x=(1/sqrt(K))*fft(y,K,1);

