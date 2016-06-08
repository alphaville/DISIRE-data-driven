function P = perm_mat(p)
%PERM_MAT returns a permutation matrix given a permutation vector
%
%P = perm_mat(p)
%
%P is such that
%
% A(:, p) = A*P
%
%and also
%
% A(p, :) = P*A
%
%It is also interesting to notice that if p and q are permutation vectors
%and
%
% P = perm_mat(p);
% Q = perm_mat(q);
%
%then it holds that
%
% P*Q = perm_mat(q(p));
%

narginchk(1,1);
n = length(p);
P = zeros(n,n);
for i=1:n,
    P(i,p(i))=1;
end