function res = AA_P_BPDN_DOA(y, A, B, r, epsilon, maxiter, tol)

% res = AA_P_BPDN(y, A, B, r, epsilon, maxiter, tol)
% 
% AA_P_BPDN, alternating algorithm for P-BPDN, solves the perturbed BPDN
% problem
%    min ||x||_1, 
%    subject to ||y - (A + B * diag(beta)) * x||_2 < epsilon
%    with respect to x, beta \in [-r,r]^n
% using an alternating approach
% 
% Written by Zai Yang, 18 Sept., 2011
% 
% References: 
% Z. Yang, C. Zhang, and L. Xie, "Stable signal recovery 
% in compressed sensing with a structured matrix perturbation", 
% 2012 IEEE International Conference on Acoustics, Speech and 
% Signal Processing (ICASSP), pp. 2737--2740, March 2012.
% 
% Z. Yang, C. Zhang, and L. Xie, "Robustly stable signal recovery 
% in compressed sensing with structured matrix perturbation", 
% IEEE Transactions on Signal Processing, vol. 60, no. 9, pp. 4658--4671, 2012.

if nargin < 6
    maxiter = 200;
end
if nargin < 7
    tol = 1e-6;
end

[m,n] = size(A);

beta = zeros(n,1);
x = zeros(n,1);

seq_x_1norm = zeros(maxiter,1);

for iter = 1:maxiter
    x_last = x;
    
    Phi = A + B * diag(beta);
    
    % solve x with beta being fixed
    cvx_begin
       variable x(n) complex;
       minimize ( norm(x,1) );
       subject to 
          norm(y - Phi * x) <= epsilon;
    cvx_end
    
    seq_x_1norm(iter) = norm(x,1);
    
    % solve beta with x being fixed
    Bdiagx = B * diag (x);
    resid = y - A * x;
    cvx_begin
       variable beta(n);
       minimize ( norm(Bdiagx * beta - resid) );
       subject to
          beta >= -r * ones(n,1);
          beta <= r * ones(n,1);
    cvx_end
    
    if abs(norm(x,1) - norm(x_last,1)) / (norm(x_last,1) + eps) < tol
        break
    end
end

res.x = x;
res.beta = beta;
res.iter = iter;

res.seq_x_1norm = seq_x_1norm(1:iter);

end