function x = sparse_sensing12( A,y,epsE )
% x = sparse_sensing12( A,y,epsE )
%
% Written by Prof. Yoash Levron, 
% Electrical Engineering, Technion, Israel, September 2014.
%
% This function solves the underdetermined system of equations Ax=y,
% with a matrix A that has less rows than columns. The function
% locates the solution vector x which is "the most sparse", that is,
% the solution vector that has the fewest number of non-zero elements. 
%
% function inputs:
% A - the sensing matrix  (dimensions M x N, where M<N)
% y - the known output vector (vector of known measurements) (dimensions M x 1)
% epsE - tolerable error of solution. If the solution vector x generates 
% an error with a second norm larger than epsE, the function throws an
% error message.
%
% function output:
% x - the sparse vector that is estimated (dimensions N x 1).
% It is assumed that x has at most 2 non-zero elements.
%
% This version of the code solves for x only if it has 1 or 2 non-zero elements,
% and it does so by scanning all possible support bases for the vector x.
% (that is, it scans every possible combination of one or two column vectors
% from A). If there is no sparse solution x (with one or two elements),
% the function will generate an error, Although a solution vector x with
% more than two non-zero elements may exist. (which can be easily checked
% by typing x = A\y in the Matlab command window)
% Because of this exhaustive scan, if a solution exists, the function
% is guranteed to locate it, and further more, the function will locate
% the sparse solution for which the second norm of the error is minimal.
% In this regard this function is much more succesful than other
% sparse solution methods such as Orthogonal Matching Pursuit, which
% are not guranteed to converge in the general case. For example, the
% matlab '\' operator (A\y) fails in the general case to find the
% most sparse solution.
%
% Conditions on the matrix A:
% The function will locate a sparse solution, if it exists. However, 
% in general, this solution may not be unique, that is, other sparse
% solutions to the system of equations Ax=y may exist.
% A sufficient condition for the solution to be unique is that every
% combination of M columns from the matrix A is linearly independant,
% where M is the number of rows in A. (The general theorem is given by
% a condition on the 'spark' of the matrix A)
%
% Complexity limit:
% The complexity of computation is determined primarily by the number
% of columns in A. The function can easily handle a matrix with
% 500 columns on a modern PC computer.


[M N] = size(A);

% check the input
there_is_an_error = 0;
if (~(size(y,1) == M))
    there_is_an_error = 1;
elseif (~(size(y,2)==1))
    there_is_an_error = 1;
elseif (M<2)
    there_is_an_error = 1;
elseif (M>N)
    there_is_an_error = 1;
end

if (there_is_an_error)
    disp('error in sparse_sensing12: mismatch in input variables');
    beep;  x = NaN;   return
end

% for the unprobable case of a square matrix :
if (M==N)
    detA = det(A);
    if ( (abs(detA))^(1/M) > epsE )
        % there is only one solution ...
        x = A\y;
        return
    end
end

% compute the second norm of each column vector
col_norm2=(sum(A.^2,1)).';   % squared second norm
col_norm = col_norm2.^0.5; % second norm

% create D, a normalized dictionary,
% which contains normalized column vectors from A
TT = repmat((col_norm.'),M,1);
D = A./TT;

% try to locate a solution vector x with a single non-zero element.
% computes the inner products with the output vector:
cc = (D.')*y;

% Choose the column for which the projection of the output is maximal:
x = zeros(N,1);
[cc_max, ii] = max(abs(cc));
ii=ii(1);
x(ii) = cc(ii);

err_vec = y-x(ii)*D(:,ii);
err = ( (err_vec.')*err_vec )^0.5;

if (err>epsE)
    % In case a sparse solution with one non-zero element does not exist,
    % try to locate a solution vector x with two non-zero elements.
    x = zeros(N,1);
    
    % scan every combination of two column vectors from A,
    % and generate an error matrix.
    err_mat = inf*ones(N,N);
    for jj1 = 1:(N-1),
        for jj2 = (jj1+1):N,
            ua = D(:,jj1);
            ub = D(:,jj2);
            U = [ua, ub];
            
            % find a minimum-energy solution to the linear system U*s = y.
            % if the columns are linearly dependant the function issues
            % a warning concerning the matrix rank. In this case the
            % function will locate a sparse solution, if it exists, but
            % this solution may not be unique.
            s = U\y;
            err_vec = U*s-y;
            err_mat(jj1,jj2) = (err_vec.')*err_vec;
        end
    end
    % locate the base which produces minimal error:
    [minerr,ind] = min(err_mat(:));  ind=ind(1);
    [opt_jj1 opt_jj2]= ind2sub(size(err_mat), ind);
    % solve again with this base to recover x:
    ua = D(:,opt_jj1);
    ub = D(:,opt_jj2);
    U = [ua, ub];
    s = U\y;
    err_vec = U*s-y;
    err = ( (err_vec.')*err_vec )^0.5;
       
    x(opt_jj1) = s(1);
    x(opt_jj2) = s(2);
         
    if (err>epsE)
        there_is_an_error = 1;
    end
end

if (there_is_an_error)
    disp('error in sparse_sensing12: a sparse solution was not found');
    x = NaN;   return
else
    % "stretch" x to match the the columns in the original matrix A
    x = x./col_norm;
end

end



