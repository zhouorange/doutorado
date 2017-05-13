function [x,r,residual, error] = omp( A, b, k,opts)
% x = OMP( A, b, k )
%   uses the Orthogonal Matching Pursuit algorithm (OMP)
%   to estimate the solution to the equation
%       b = A*x     (or b = A*x + noise )
%   where there is prior information that x is sparse.
%
%   "A" may be a matrix, or it may be a cell array {Af,At}
%   where Af and At are function handles that compute the forward and transpose
%   multiplies, respectively.
%
% [x,r,residual,error] = OMP( A, b, k, opt )
%   is the full version.
% Outputs:
%   'x' is the k-sparse estimate of the unknown signal
%   'r' is the residual b - A*x
%   'normR' = norm(r)
%   'residual'     is a vector with normR from every iteration
%   'error'       is a vector with the output of errFcn from every iteration
%
% Inputs:
%   'A'     is the measurement matrix
%   'b'     is the vector of observations
%   'k'     is the sparsity
%
%   'opts'  is a structure with more options, including:
%       .target     = Define a residual in case you do not have knowladge
%                     of the sparsity.
%       .mode   = {0,1};  
%       If mode=0, there is no orthogonalization procedure. 
%       If mode=1 it is performed an orthonormalization on the matrix A.
%

if nargin<4
    target_resid    = 1e-12;     
    mode            = 0    ;
elseif nargin == 4
     mode = opts.mode;
     target_resid = opts.target;   
end

% 
if iscell(A)
    LARGESCALE  = true;
    Af  = A{1};
    At  = A{2};    
else
    LARGESCALE  = false;
    Af  = @(x) A*x;
    At  = @(x) A'*x;
end


%%                                        -- Intitialize --
% start at x = 0, so r = b - A*x = b
r           = b;
normR       = norm(r);
Ar          = At(r);
N           = size(Ar,1);       
M           = size(r,1);        
if k > M
    error('K cannot be larger than the length sequence');
end
unitVector  = zeros(N,1);
x           = zeros(N,1);

indx_set    = zeros(k,1);
indx_set_sorted     = zeros(k,1);
A_T         = zeros(M,k);
A_T_nonorth = zeros(M,k);
residual   = zeros(k,1);
error     = zeros(k,1);



for kk = 1:k
    
    [~,ind_new]     = max(abs(Ar));
    
    indx_set(kk)    = ind_new;
    indx_set_sorted(1:kk)   = sort( indx_set(1:kk) );
    
    if LARGESCALE
        unitVector(ind_new)     = 1;    [~,index(ii)] = max(abs(Ar'*r));

        atom_new                = Af( unitVector );
        unitVector(ind_new)     = 0;
    else
        atom_new    = A(:,ind_new);
    end
    
    A_T_nonorth(:,kk)   = atom_new;       
    
    
    % -- Step 2: update residual
    
    if mode==0
        % The straightforward way:
        x_T = A_T_nonorth(:,1:kk)\b;    
      
        x( indx_set(1:kk) )   = x_T;
        r   = b - A_T_nonorth(:,1:kk)*x_T;
    elseif mode==1
    
        % We use Gram-Shimidt
        for j = 1:(kk-1)
            atom_new    = atom_new - (A_T(:,j)'*atom_new)*A_T(:,j);
        end
        % Second, normalize:
        atom_new        = atom_new/norm(atom_new);
        A_T(:,kk)       = atom_new;
        % solve least-squares problem 
        x_T     = A_T(:,1:kk)'*b;
        x( indx_set(1:kk) )   = x_T;      % note: indx_set is guaranteed to never shrink
        % Fourth, update residual:
        r       = b - A_T(:,1:kk)*x_T;

    end
    
    
    normR   = norm(r);

    residual(kk)   = normR;
    

    
    if normR < target_resid
        if PRINT
            fprintf('Residual reached desired size (%.2e < %.2e)\n', normR, target_resid );
        end
        break;
    end
    
    if kk < k
        Ar  = At(r); 
    end
    
end

if mode==1 
 % For the last iteration, we need to do this without orthogonalizing A
 % so that the x coefficients match what is expected.
 x_T = A_T_nonorth(:,1:kk)\b;
 x( indx_set(1:kk) )   = x_T;
end
r       = b - A_T_nonorth(1:kk)*x_T;
normR   = norm(r);

end % end of main function