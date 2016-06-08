function Fg = F_gamma(A,x,y,tau,d,gamma,lambda)

xtd = x+tau*d;

nabla = A'*(A*x-y);
nabla_xtd = A'*(A*xtd-y);

Fg_1 = 0.5*(A*xtd-y)'*(A*xtd-y) - 0.5*gamma*(nabla_xtd'*nabla_xtd) + lambda*g_gamma(xtd-gamma*nabla_xtd, gamma*lambda);
Fg_2 = 0.5*(A*x-y)'*(A*x-y)     - 0.5*gamma*(nabla'*nabla)         + lambda*g_gamma(x-gamma*nabla, gamma*lambda);
Fg = Fg_1 - Fg_2;

