function gg = g_gamma(x, gl)
gg = 0;
for i=1:length(x),
    xi = abs(x(i));
    if xi <= gl,
        gg = gg + xi^2/(2*gl);
    else
        gg = gg + xi - (gl/2);
    end
end
