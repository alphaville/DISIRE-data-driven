function e = l1_menv(x, gamma)
%L1_MENV is the Moreau envelope of l1-norm.


n = length(x);
e = 0.0;
for i=1:n,
    xi = abs(x(i));
    if xi <= gamma
        e = e+ xi^2/(2*gamma);
    else
        e = e + xi - (gamma/2);
    end
end