function [Y,u,freq,amp] = ANM_sdpt3(YonOmega, Omega, N, eta)

% [Y,u,freq,amp] = ANM_sdpt3(YonOmega, Omega, N, eta)
% 
% ANM_sdpt3 implements the atomic norm minimization problem:
% min_Y ||Y||_A, subject to ||Y_Omega - Yo_Omega||_F <= epsilon
% via duality using SDPT3.
% 
% The dual problem is
% min_{V,H} <V_Omega, Yo_Omega>_R + epsilon*||V_Omega||_F,
% subject to [I V'; V H] >= 0, V_Omegac = 0, and T^*(H) = [1,0,...,0]^T.
% 
% Input:
%   YonOmega: observed measurements on Omega
%   Omega: index set of the measurements
%   N: length of sinusoidal signal of interest
%   eta: upperbound of the Frobenius norm of noise
% Output:
%   Y: recovered Y
%   u: u composing the Toeplitz matrix
%   freq: recovered frequency
%   amp: recovered amplitude
% 
% References:
% Z. Yang and L. Xie, "Continuous Compressed Sensing With a Single or ...
%     Multiple Measurement Vectors", IEEE Workshop on Statistical Signal ...
%     Processing (SSP), pp. 308--311, June 2014.
% Z. Yang and L. Xie, "Exact joint sparse frequency recovery via ...
%     optimization methods", http://arxiv.org/abs/1405.6585, May 2014.
% 
% Written by Zai Yang, Feb. 2014


if nargin < 4
    eta = 0;
end
if nargin < 3
    N = Omega(end);
end

L = size(YonOmega,2);
Omegac = (1:N)';
Omegac(Omega) = [];


% solve the dual problem

cvx_quiet false
cvx_precision default
cvx_solver sdpt3
cvx_begin sdp 

  variable X(N+L,N+L) hermitian,
  dual variable U,
  
  X >= 0 : U,
  X(1:L,1:L) == eye(L),
  X(L+Omegac,1:L) == 0,
  trace(X) == L + 1,
  for j = 1:N-1
      sum(diag(X(L+1:N+L,L+1:N+L), j)) == 0;
  end
  
  if eta == 0
      minimize real(trace(X(L+Omega,1:L)'*YonOmega));
  else
      minimize real(trace(X(L+Omega,1:L)'*YonOmega)) + ...
          eta*norm(X(L+Omega,1:L),'fro');
  end
  
cvx_end

% Y estimate
Y = U(L+1:L+N,1:L) * 2;

% postprocessing
u = U(L+1,L+1:L+N).';
[freq, amp] = VanDec(u);

amp = amp * 2/sqrt(L);

end