function A = ambiguity_function(gamma,g,tau_max,Td)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%     Generates Ambiguity Function


G_matr=conj(toeplitz([g zeros(1,Td-length(g))],[g(1) zeros(1,tau_max-1)]));
Gamma_G_matr=repmat([gamma zeros(1,Td-length(gamma))].',1,tau_max).*G_matr;

A=fft(Gamma_G_matr);