function [h, g] = ldpc_h2g(H)
% converts tentative binary LDPC matrix H into a new matrix h
% (columns are permuted) and produces the generator matrix g
% H should be a sparse matrix in MATLAB format.
%
% MEX file


%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.0 $  $Date: 1999/08/23 $
