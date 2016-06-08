function [x_star, details] = lasso_newton2(A, y, lambda, x0, ops)
%LASSO_NEWTON solves the LASSO problem
%
%  min_x ||Ax - y||^2 + lambda*||x||_1
%
%

% Debug mode
dbg = true;



% Processing input
narginchk(3,5);
default_ops = struct('max_iter', 1000, 'tol', 1e-6, 'sigma', 1e-4);
if nargin<=4,
    ops = default_ops;
    if nargin==3,
        x0 = zeros(size(A,2),1);
    end
else
    if ~isfield(ops,'max_iter'),
        ops.max_iter = default_ops.max_iter;
    end
    if ~isfield(ops,'tol'),
        ops.tol = default_ops.tol;
    end
    if ~isfield(ops,'sigma'),
        ops.sigma = default_ops.sigma;
    end
end
max_iter = ops.max_iter;
tolerance = ops.tol;
sigma = ops.sigma;




% Constants
tic;
q=(y'*A)';
[m,n] = size(A);
Hess = @(x) A'*(A*x);
eigsOpt.issym = 1;
eigsOpt.tol = 1e-3;
Lf = eigs(Hess, n, 1, 'LM', eigsOpt);
gamma = 0.95/Lf;




% Initialization
x = x0;
res = A*x - y;          % res:  residual (Ax-y)
zeta = (res'*A)';       % zeta: gradient = A'*res = A'(Ax-y)
L=[];
alpha_p = [];           % alpha_p : previous set of alpha
p = min(size(A));       % p: smallest dimension of A (rows)
if dbg,
    cost = zeros(max_iter,1);
    taus = zeros(max_iter,1);
end
d_alpha = [];           % set of rejected indices from alpha
lambda_target = lambda;

for i = 1:max_iter,
    if dbg, cost(i) = 0.5*norm(res,2)^2 + lambda*norm(x,1); end
    h = x - gamma*zeta;    
    absh_over_gam = abs(h)/gamma;
    alpha = find(absh_over_gam > lambda)';
    alpha_c = setdiff(1:n, [alpha d_alpha]);
    % compute residual
    S = sign(h);
    r = zeros(n,1);
    r(alpha)   = gamma*(zeta(alpha)   + lambda * S(alpha)  );
    r(d_alpha) = gamma*(zeta(d_alpha) + lambda * S(d_alpha) );
    r(alpha_c) = x(alpha_c);
    resnorm = norm(r/gamma);
    if dbg, resnorms(i) = resnorm; end
%     if (lambda == lambda_target) && resnorm < 1e-8
%         break
%     elseif resnorm < 1e-8
%         lambda = 0.5*lambda;
%     end
% 
%     if length(alpha) > m
%         [hsort,isort] = sort(absh_over_gam,1,'descend');
%         lambda = hsort(m+1)/gamma;
%         alpha = sort(isort(1:m))';
%         alpha_c = setdiff(1:n, [alpha d_alpha]);
%     end
    
    % alpha is the set of indices i where |x(i) - gamma * gradient_i f(x)|
    % = |x(i) - gamma*zeta(i)| > gamma*lambda
    %     alpha = find(abs(h) > gamma*lambda)';
    
%     % Construct of update L
%     if isempty(L),
%         % When L is not available (first iteration)...
%                 
%         % The first time, we need to use `chol` to factorize A. We
%         % therefore restrict the size of alpha to p (smallest size of A (rows))
%         a_len = length(alpha);
%         if (a_len > p),
%             d_alpha = alpha(p+1:end);
%             alpha=alpha(1:p);
%         end
%         L = chol(A(:, alpha)'*A(:, alpha), 'lower');
%     else
%         alpha_before = alpha;
%         [L, alpha] = chol_ata_update(A, L, alpha_p, alpha_before);
%         d_alpha = setdiff(alpha_before, alpha);        
%     end
%     
%     
%     if dbg,
%         %Verify the Cholesky factorization and check the diagonal entries
%         %of L (they must be positive)
%         if~(norm(L*L' - A(:,alpha)'*A(:,alpha), Inf)<1e-12),keyboard,end
%         if~(min(diag(L)) > 1e-6),keyboard,end              
%     end
%     
%     
%     % Here, alpha is the set original set minus those indices that make
%     % A(:, alpha)'*A(:, alpha) signular. These are in d_alpha.
%     
%     % We define alpha_c as all other indices (not in alpha and d_alpha)
%     
    if dbg,
        %Also make sure that alpha, d_alpha and alpha_c make sense:
        assert( length(unique([alpha alpha_c d_alpha])) == n );
    end
    
        
    % r is R_gamma(x)
    %

    
  % Compute direction d
    d = zeros(length(x),1);
    d(alpha_c) = -x(alpha_c);
    d(d_alpha) = -r(d_alpha);
    b = -(zeta(alpha)+lambda*sign(h(alpha)));
%     if ~isempty(d_alpha),
%         b = b - A(:, alpha)'*A(:, d_alpha)*d(d_alpha) ...
%             - A(:, alpha)'*A(:, alpha_c) * d(alpha_c);
%     end
%     d(alpha) = L'\(L\b);
    d(alpha)  = (A(:, alpha)'*A(:, alpha)+1e-6*resnorm*eye(length(alpha)))\b-x(alpha);
%     d(alpha)  = (A(:, alpha)'*A(:, alpha))\b-x(alpha);
    
%     if dbg,
%         assert( norm (d(alpha)   -   L'\(L\b), Inf) < 1e-6);
%     end
    
    
    % Line search (Armijo):
    xi = A*d;
    delta = (xi'*A)';
    s = d - gamma*delta;
    beta_1 = s'*zeta;
    beta_2 = 0.5*(s'*delta);
    rho = (sigma/gamma)* (r'*s);
    k = lambda*g_gamma(h, gamma*lambda);
    tau = 1;
    epsilon = 0.5;
    do_run = true;
    if rho>0
        keyboard
    end
    while do_run && tau > 1e-12,
        tau = epsilon*tau;
        do_run = (beta_1 * tau + beta_2 * tau^2 + ...
            lambda*g_gamma(h + tau * s, gamma*lambda) - k > tau*rho);
    end
    
    % Updates of zeta, x and res
    if dbg, taus(i) = tau; end
    zeta = zeta + tau * delta;
    x    = x    + tau * d;
    res  = res  + tau * xi;
    
    %if (norm(tau*d, Inf) < tolerance), break; end
    
    alpha_p = alpha;
end
x_star = x;


% Output some details
details.iter=i;
if dbg,  details.cost = cost(1:i,1);  end
if dbg,  details.resnorm = resnorms;  end
details.flag=0;
details.time = toc;
details.tau = taus;
details.tolerance = tolerance;

