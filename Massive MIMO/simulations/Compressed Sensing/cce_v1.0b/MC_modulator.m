function s = MC_modulator(a,g,N)

%%%    Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%	Multicarrier Modulator (Pulse-Shaping)
%%%
%%%%%%%%%%%%%%%%%%%
%%%
%%%   a   ...     matrix consisting of data symbols 
%%%               number of rows = number of sub-carriers (K)
%%%               number of columns = number of transmitted data symbols (L)
%%%   g   ...     transmit pulse (vector)
%%%               assumed to be causal, i.e., 
%%%               g=(g[0] ... g[length_g-1]), otherwise zero
%%%               in MATLAB - of course - g=g(1:length_g) , i.e., a row
%%%               vector
%%%   N   ...     discrete symbol period
%%%
%%%   s   ...     transmit signal (column vector)
%%%
%%%
%%% we implement: 
%%%       s(n)=1/sqrt(K)sum_{l=0}^{L-1} sum_{k=0}^{K-1}
%%%                               a_{k,l}g(n-lN)e^{j 2 pi k/K(n-lN)}

[K L]=size(a);
length_g=length(g);

b=sqrt(K)*ifft(a, K, 1);

c=b;
for zaehler=1:(ceil(length_g/K)-1)
      c=[c;b];
end    

c=c(1:length_g,:);

g_transp=g.';

g_matr=[];
for zaehler=1:L
    g_matr=[g_matr g_transp];
end    

d=g_matr.*c;

no_of_rows_desired=ceil(length_g/N)*N;
d=[d ; zeros((no_of_rows_desired-length_g),L)];

no_of_ISI_symbs=ceil(length_g/N)-1;

d=[zeros(no_of_rows_desired,no_of_ISI_symbs) d zeros(no_of_rows_desired,no_of_ISI_symbs)];

e=[];
for zaehler=1:(L+no_of_ISI_symbs)
    spalt_vec=zeros(N,1);
    for zaehler1=0:no_of_ISI_symbs
        spalt_vec=spalt_vec+d((1+N*zaehler1):(N+N*zaehler1),zaehler+no_of_ISI_symbs-zaehler1);
    end    
    e=[e spalt_vec];
end    

s = e(:);
