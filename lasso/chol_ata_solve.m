function x = chol_ata_solve(L, b, ptr)
%CHOL_ATA_SOLVE solves the Cholesky system LL'x = b or the perturbed
%Cholesky system PLL'P'x = b.
%
%
if nargin==2,
    x = (L')\(L\b);
else
    [~, p] = sort(ptr);
    y = L'\(L\b(ptr));
    x = y(p);
end