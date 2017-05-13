% This code demonstrates the use of the sparse sensing funtion
% function x = sparse_sensing12( A,y,epsE )
%
% (See documentation within function or type 'help sparse_sensing12'
%
% Written by Yoash Levron, 
% Electrical Engineering Faculty, Technion, Israel, September 2014.
%
% Please select an example below by commenting the other example 
% using the '%' symbol.

clc;
epsE = 1e-9;   % tolerable solution error

% %%%% example 1:  
% %%%% (recovers a sparse vector with a single non-zero element
% %%%% from a length 2 output vector.
% disp('EXAMPLE 1');
% disp('---------');
% a1 = [2 ; 0];  
% a2 = [3 ; 4];  
% a3 = [2 ; 1];
% a4 = [1 ; 2];
% a5 = [5 ; 1];
% a6 = [-3 ; 3];
% A = [a1, a2, a3, a4, a5, a6];
% 
% y = A * [0 ; 0 ; 0 ; 18 ; 0 ; 0];  % or:  y = 18*a4
% xa = sparse_sensing12( A,y,epsE );
% xa


%%% example 2:  
%%% recovers a sparse solution with a two non-zero elements
%%% from a length 4 output vector.
disp(' ');
disp('EXAMPLE 2');
disp('---------');

a1 = [2 ; -10 ; 50 ; 6];  
a2 = [3 ; 4 ; 2 ; -1];  
a3 = [2 ; 1 ; 8 ; 1];
a4 = [1 ; 2 ; 13 ; 4];
a5 = [5 ; 1 ; 9 ; 3];
a6 = [-3 ; 5 ; 2 ; 7];
a7 = [-3 ;-3 ; -2; 5];
a8 = [-4 ; -6 ; -4; 6];
a9 = [5 ; -2 ; 8 ; 8];
a10 = [4 ; 5 ; -1 ; -2];
A = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10];

y = 12*a4 + (-7)*a9;
xb = sparse_sensing12( A,y,epsE );
xb


% %%% example 3:  
% %%% recovers a sparse vector with two non-zero elements
% %%% from a length 4 output vector.
% %%% The matrix A is now randomized.
% disp(' ');
% disp('EXAMPLE 3');
% disp('---------');
% 
% A =20*rand(4,14) - 10;
% 
% a6 = A(:,6);   % 6'st column
% a11 = A(:,11);   % 11'st column
% 
% y = 5*a6 + 9*a11;
% xc = sparse_sensing12( A,y,epsE );
% xc


% %%% example 4:  
% %%% the function fails to locate a sparse solution with three
% %%% non-zero elements (present version of the function only handles two
% %%% non-zero elements.
% disp(' ');
% disp('EXAMPLE 4');
% disp('---------');
% 
% A =20*rand(6,14) - 10;
% 
% a3 = A(:,3);   % 3'rd column
% a6 = A(:,6);   % 6'st column
% a11 = A(:,11);   % 11'st column
% 
% y = (-6)*a3 + 5*a6 + 9*a11;
% xd = sparse_sensing12( A,y,epsE );
% xd




