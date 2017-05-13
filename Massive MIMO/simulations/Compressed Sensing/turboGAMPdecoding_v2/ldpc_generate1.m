function [H] = generate(M,N,t,q,seed)
%[H] = generate(M,N,t,q,seed)
% generates sparse LDPC parity check matrix over GF2
% Input: 
%  M      number of parity checks
%  N      blocklength
%  t      mean column weight
%  q      GF base (only 2 now)
%  seed   initializes random generator - we use MATLAB 5 uniform generator
%MEX file
% C code is based on Davey&MacKay code http://wol.ra.phy.cam.ac.uk/mackay/codes/
%   (sparse.c) with the author's permission
%
% Example
%           [H] = generate(4512,6112,2.3,2,123);

%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
