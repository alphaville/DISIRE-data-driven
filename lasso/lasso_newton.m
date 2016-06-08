function [x_star, details] = lasso_newton(A, y, lambda, x0, ops)
%LASSO_NEWTON solves the LASSO problem
%
%  min_x ||Ax - y||^2 + lambda*||x||_1
%
%

% Debug mode
dbg = true;



% Processing input
narginchk(3,5);
default_ops = struct('max_iter', 100, 'tol', 1e-6, 'sigma', 1e-4);
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
gamma = 0.997/(normest(A,2)^2);
n = size(A,2);




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



for i = 1:max_iter,
    if dbg, cost(i) = 0.5*norm(res,2)^2 + lambda*norm(x,1); end
    h = x - gamma*zeta;
    
    % alpha is the set of indices i where |x(i) - gamma * gradient_i f(x)|
    % = |x(i) - gamma*zeta(i)| > gamma*lambda
    alpha = find(abs(h) > gamma*lambda)';
    
    % Construct of update L
    if isempty(L),
        % When L is not available (first iteration)...
                
        % The first time, we need to use `chol` to factorize A. We
        % therefore restrict the size of alpha to p (smallest size of A (rows))
        a_len = length(alpha);
        if (a_len > p),
            d_alpha = alpha(p+1:end);
            alpha=alpha(1:p);
        end
        L = chol(A(:, alpha)'*A(:, alpha), 'lower');
    else
        alpha_target = alpha;
        [L, alpha] = chol_ata_update(A, L, alpha_p, alpha_target);
        d_alpha = setdiff(alpha_target, alpha);        
    end
    
    
    rcond(A(:, alpha)'*A(:, alpha) + 1e-8*eye(length(alpha)))
    
    if dbg,
        %Verify the Cholesky factorization and check the diagonal entries
        %of L (they must be positive)
        assert(norm(L*L' - A(:,alpha)'*A(:,alpha), Inf)<1e-12);
        assert(min(diag(L)) > 1e-6)                
    end
    
    
    % Here, alpha is the set original set minus those indices that make
    % A(:, alpha)'*A(:, alpha) signular. These are in d_alpha.
    
    % We define alpha_c as all other indices (not in alpha and d_alpha)
    alpha_c = setdiff(1:n, [alpha d_alpha]);
    
    if dbg,
        %Also make sure that alpha, d_alpha and alpha_c make sense:
        assert( length(unique([alpha alpha_c d_alpha])) == n );
    end
    
    S = lambda * sign(h(alpha));
    
    
    % r is R_gamma(x)
    %
    r = zeros(n,1);
    r(alpha)   = gamma*(zeta(alpha)   + S  );
    r(d_alpha) = gamma*(zeta(d_alpha) + lambda * sign(h(d_alpha)) );
    r(alpha_c) = x(alpha_c);    
    
    
  % Compute direction d
    d = zeros(length(x),1);
    d(alpha_c) = -x(alpha_c);
    d(d_alpha) = -gamma*(zeta(d_alpha) + lambda * sign(h(d_alpha)));
    b = -(zeta(alpha)+lambda*sign(h(alpha)));
    if ~isempty(d_alpha),
        b = b - A(:, alpha)'*A(:, d_alpha)*d(d_alpha) ...
            - A(:, alpha)'*A(:, alpha_c) * d(alpha_c);
    end
    %d(alpha) = L'\(L\b);
    d(alpha)  = (A(:, alpha)'*A(:, alpha))\b;   
    
    
    % Line search (Armijo):
    xi = A*d;
    delta = (xi'*A)';
    s = d - gamma*delta;
    beta_1 = s'*zeta;
    beta_2 = 0.5*(s'*delta);
    rho = (sigma/gamma)* (r'*s);
    k = lambda*g_gamma(h, gamma*lambda);
    tau = 2;
    epsilon = 0.5;
    do_run = true;   
    
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
details.flag=0;
details.time = toc;
details.tau = taus;
details.tolerance = tolerance;

