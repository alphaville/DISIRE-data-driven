function [l1, l2, P_tr, flag]=chol_ata_insert_col(L, A, c, idx)
%CHOL_ATA_INSERT_COL is similar to chol_ata_append_col, but the difference
%is that here the additional column is inserted in A. Let A be a matrix
%in the following form
%
% A = [a(1) ... a(idx-1) a(idx) ... a(n)]
%
%where a(i) are column-vectors. We insert a column-vector c to get the
%matrix
%
% A_ = [a(1) ... a(idx-1) c a(idx) ... a(n)]
%
%This function computes the permuted Cholesky factorisation of A_'*A_ given
%the Cholesky factorisation of A'*A, that is, it constructs a lower
%triangular matrix L_ and a permutation matrix P_ so that
%
% A_'A_ = P_ L_ L_' P_'
%
%The matrix L_ has the following form
%
% L_ = [ L    0
%        l1'  l2]
%
%The (transpose of the) permutation matrix P_' is returned as a vector.
%Notice that P_' is such that A_*P_' = [A c]. Using the returned vector of
%integer indices P_, the following holds true:
%
% Ac'*Ac == T(P_tr, P_tr))
% 
%Syntax:
% [l1, l2, P_tr, flag]=CHOL_ATA_INSERT_COL(L, A, c, idx)
%
%Input arguments:
% L         Lower triangular matrix. The Cholesky factorisation of A'*A.
% A         Original matrix A
% c         Column to add in A
% idx       Index where c will be added
%
%Output arguments:
% l1, l2    Data used to update the Cholesky factorisation of A'A
% P_tr      Transpose of the permutation matrix P as described above. Note
%           that it is more convenient to return P_' (= P_tr) than P. To 
%           retrieve the permutation vector P_, if necessary, use the
%           command [~, P_] = sort(P_tr).
% flag      this flag is set to 0 if the computation has succeded and to 
%          -1 if the new matrix [A c]'[A c] is not positive definite.
%
%See also:
% chol_ata_append_col, chol_ata_remove_col, chol_ata_solve

% Pantelis Sopasakis

% GIST:
% https://gist.github.com/alphaville/4459f416c3d790b43502

narginchk(4,4);
n = size(A, 2);
P_tr = [1:idx-1 n+1 idx:n];
[l1, l2, flag] = chol_ata_append_col(L,A,c);