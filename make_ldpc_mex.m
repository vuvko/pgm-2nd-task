% MAKE_LDPC_MEX is a function to generate a parity check matrix for LDPC code
%
% Matlab interface to Software for Low Density Parity Check Codes by
% Radford M. Neal and Peter Junteng Liu:
% http://www.cs.utoronto.ca/~radford/ftp/LDPC-2012-02-11/index.html
%
% Usage:
% matrix = make_ldpc_mex(M, N, K);
% 
% Inputs:
% M	- number of checks of in the code (a single positive number)
% N - length of a transmitted codewords (a single positive number)
% K - desirable number of "ones" in each column of the parity check matrix (a single number)
%
% Outputs:
% matrix - the generated parity check matrix (sparse matrix of size M*N)
% 
% To build the code in matlab choose reasonable compiler ("mex -setup") and run build_make_ldpc_mex.m
% File "randfile" is required by the library at the runtime
% example:
%       checkMatrix = make_ldpc_mex(12, 16, 3);
%       disp(full(checkMatrix));
% 
% by Anton Osokin, firstname.lastname@gmail.com, Feb 2013