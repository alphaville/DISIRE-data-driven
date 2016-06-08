function result = FBNewtonL1LS(A, b, lambda0, opt)
%FBNEWTONL1LS solves a l1-regularized least squares problem with the Forward-Backward
%Newton method presented in [1]. It solves the following problem:
%
% minimize_x (1/2) ||Ax-b||^2 + lambda0 * ||x||_1
%
%
%Syntax:
% result = FBNewtonL1LS(A, b, lambda0, opt)
%
%Input arguments:
% A, b:         problem data (matrix A, vector b)
% lambda0:      regularization parameter
% opt:          options
%
% opt is a structure containing the method options; these are:
%  - direction   : 1 (Newton with LDL solves)
%                : 2 (L-BFGS)
%                : 3 (Newton with CG solves)
%
%  - linesearch  : 1 (Armijo backtracking)
%                : 2 (Wolfe conditions via bisection)
%
% It is recommended that direction 1 or 3 is used (use 3 when away from 
% the solution, and 1 when you are close to the solution, although, both
% work fine) and linesearch 1 is recommended.
%
%  - maxit       : maximum number of iterations
%                  default value: 1000
%
%  - tol         : maximum solution tolerance 
%                  default value: 1e-8
%
%  - tolc        : tolerance used for the continuation algorithm
%                  continuation happens when norm(FPR) <= lambda*opt.tolc
%                  where FPR is the fixed-point residual
%                  default value: 1e-4
%
%  - cont        : (true/false) whether to use continuation
%                  It is recommended to use the continuation method
%                  default value: true
%
%  - eta         : continuation parameter
%                  default value: 0.5
%
%  - memory      : when the direction is 2 (L-BFGS), `memory` specifies the buffer length
%                  default value: 5
%
%  - chol_update : (true/false) whether to update the Cholesky factorization
%                  default value: true
%
%  - gamma       : proximal operator parameter gamma
%                  if not provided, the solver assigns to it the value 0.95/Lf, where
%                  Lf is the gradient Lipshitz constant of f(x) = ||Ax-b||^2
%
%  - x           : initial guess (x0)
% 
%  - res         : residual of x (res = Ax - b)
%
%  - loss        : cost of the initial guess, that is ||Ax-b||^2 + lambda0 * ||x||_1
%
%  - gradf       : gradient of f at the initial guess, that is A'*res
%
%  - gs          : value of x-gamma*gradf
%  
%
%
%Output arguments:
% result         : structure which packs together all the necessary information
%                  which is
%
%  - x           : optimizer
%
%  - iterations  : number of iterations
%
%  - pobj        : the sequence of the value of the primal objective, that is 
%                  ||Ax-b||^2 + lambda0 * ||x||_1.
%                  You can plot pobj to gain some insight about how (fast) the algorithm
%                  converges
%
%  - times       : an array of runtimes per iteration
%
%  - tau         : value of tau (linesearch)
%
%  - dim         : number of nonzeros in alpha (active indices)
%
%  - cgIters     : when a conjugate gradient method is used internally, `cgIters` reports
%                  the number of iterations needed for convergence
%
%  - res, 
%  - loss, 
%  - gradf, 
%  - gamma, 
%  - gs          : data which are provided by the algorithm as output and can be provided 
%                  to the next instance of the solver
%
%  - L           : Lipschitz constant of the gradient of f(x) = ||Ax-b||^2
%
%  -alpha,beta   : active/inactie indices
%
%
%
% [1] P. Sopasakis, N.Freris and P. Patrinos (2016), Accelerated reconstruction of a 
%     compressively sampled data stream, 24th European Signal Processing conference, submitted.
%


[m,n] = size(A);
At = A';
Atb =At*b;
if nargin < 4 || ~isfield(opt, 'gamma') || isempty(opt.gamma)
    Hess = @(x) A'*(A*x);eigsOpt.issym = 1;eigsOpt.tol = 1e-5;
    Lf = eigs(Hess, n, 1, 'LM', eigsOpt);
    gamma = 0.95/Lf;
else
    gamma = opt.gamma;
end

%% ******************************************* CHECK ARGUMENTS
if nargin < 4 || ~isfield(opt, 'maxit') || isempty(opt.maxit),opt.maxit = 1e3; end
if nargin < 4 || ~isfield(opt, 'tol')   || isempty(opt.tol),opt.tol = 1e-8;  end
if nargin < 4 || ~isfield(opt, 'tolc')  || isempty(opt.tolc),opt.tolc = 1e-4;end
if nargin < 4 || ~isfield(opt, 'cont')  || isempty(opt.cont),opt.cont = true;end
if nargin < 4 || ~isfield(opt, 'eta')   || isempty(opt.eta),opt.eta = 0.5;end
if nargin < 4 || ~isfield(opt, 'direction') || isempty(opt.direction),opt.direction = 1;end
if nargin < 4 || ~isfield(opt, 'memory') || isempty(opt.memory),opt.memory = 5;end
if nargin < 4 || ~isfield(opt, 'linesearch') || isempty(opt.linesearch),opt.linesearch = 1;end
if nargin < 4 || ~isfield(opt, 'chol_update') || isempty(opt.chol_update),opt.chol_update = true;end
if opt.direction == 3
    cgIters(1,opt.maxit) = 0;
    cgHow(1,opt.maxit) = 0;
    if nargin < 4 || ~isfield(opt, 'cgMaxIter') || isempty(opt.cgMaxIter),opt.cgMaxIter = n;end
end

if nargin<4 || ~isfield(opt,'x') || ~any(opt.x)
    x = sparse(zeros(n,1));
    res = -b;
    loss = 0.5*(b'*b);
    gradf = -Atb;
    gs = gamma*A'*b;
else
    x = opt.x;
    % Initial residual
    if ~isfield(opt,'res'),
        res = A*x-b;
    else
        res = opt.res;
    end
    
    % Initial loss := 0.5*||res||^2
    if ~isfield(opt,'loss'),
        loss = 0.5*(res'*res);
    else
        loss = opt.loss;
    end
    
    % Initial gradf
    if ~isfield(opt,'gradf'),
        gradf = At*res;
    else
        gradf = opt.gradf;
    end
    
    % Initial gs
    if ~isfield(opt,'gs'),
        gs = x - gamma*gradf;
    else
        gs = opt.gs;
    end
end


if opt.cont,
    if opt.direction == 1
        ls = sort(abs(gs)/gamma,1,'descend');
        lambda = max(lambda0,ls(m+1));
    else
        lambda = max(opt.eta*norm(gradf,'inf'),lambda0);
    end
else
    lambda = lambda0;
end

% Forward-Backward Step
FBS = sign(gs).*max(0,abs(gs)-lambda*gamma);
% Fixed Point Residual
FPR = (x-FBS)/gamma;
FBE = loss + lambda*norm(FBS,1) - gamma*(gradf'*FPR) + 0.5*gamma*(FPR'*FPR);
if opt.direction == 2
    gradFBE = FPR/gamma-At*(A*FPR);
    gradFBEnrm = sqrt(gradFBE'*gradFBE);
end


%% INITIALIZE STUFF
result.times = zeros(1,opt.maxit);
result.pobj = zeros(1,opt.maxit);
result.tau = zeros(1,opt.maxit);
result.dim = zeros(1,opt.maxit);
d = zeros(n, 1);
icont = 0;
flagChangedLambda = false;
start = tic();

result.its = [];

L=[];
alpha_prev = [];
alpha=[];
beta=[];

for it = 1:opt.maxit
    % compute objective wrt to target lambda
    result.pobj(1,it) = loss + lambda0*norm(x,1);
    
    % store results
    result.times(1,it) = toc(start);
    result.FPR(1,it) = norm(FPR);
    % stopping criterion
    if  lambda == lambda0 && norm(FPR) <= opt.tol
        if opt.cont
            result.its(icont+1) = it-sum(result.its(1:icont))+1;
        end
        break;
    end
    switch opt.direction
        case 1  %Newton direction
            % inactive (beta) and active (alpha) variables
            beta = find(abs(gs) <= lambda*gamma);
            alpha =find(abs(gs) > lambda*gamma);
            result.dim(1,it) = nnz(alpha);
            if nnz(alpha)>m% if continuation is used this is not needed
                reg = 1e-6*result.FPR(1,it)*speye(nnz(alpha));
                % disp('Regularization used')
            else
                reg = 0;
            end
            % solve explicitly for inactive components
            d(beta) = -x(beta);
            % solve for active components
            
            if opt.chol_update && opt.cont,
                if isempty(alpha_prev),
                    fprintf('Performing Cholesky with size(alpha) = %d\n', length(alpha));
                    cholesky_start = tic;
                    L = chol(A(:, alpha)'*A(:, alpha), 'lower');
                    fprintf('Cholesky done in %g\n', toc(cholesky_start));
                else
                    % TODO: Check whether updating the Cholesky would incur
                    % more FLOPS than doing it from scratch
                    [L, alpha] = chol_ata_update(A, L, alpha_prev, alpha);
                end
                rhs = Atb(alpha)-lambda*sign(gs(alpha));
                d(alpha) = L'\(L\rhs) - x(alpha);                
            else
                rhs = Atb(alpha)-lambda*sign(gs(alpha));
                XX = At(alpha,:)*At(alpha,:)' + reg * eye(length(alpha));
                d(alpha) = XX\rhs - x(alpha);
            end
            
            
            
            
        case 2  % L-BFGS
            if it == 1 || flagChangedLambda
                flagChangedLambda = false;
                alphaC = 0.01;
                skipCount = 0;
                d = -gradFBE; % use steepest descent direction initially
                LBFGS_col = 1;
                LBFGS_mem = 0;
            else
                %%% x' - x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Sk = x - xp;
                Yk = gradFBE - gradFBEp;
                YSk = Yk'*Sk;
                if gradFBEnrm < 1
                    alphaC = 3;
                end
                if YSk > 1e-6*gradFBEnrm^alphaC;
                    LBFGS_col = 1+mod(LBFGS_col, opt.memory);
                    LBFGS_mem = min(LBFGS_mem+1, opt.memory);
                    S(:,LBFGS_col) = Sk;
                    Y(:,LBFGS_col) = Yk;
                    YS(LBFGS_col) = YSk;
                else
                    skipCount = skipCount+1;
                end
                if LBFGS_mem > 0
                    H = YS(LBFGS_col)/(Y(:,LBFGS_col)'*Y(:,LBFGS_col));
                    d = LBFGS(S, Y, YS, H, -gradFBE, int32(LBFGS_col), int32(LBFGS_mem));
                else
                    d = -gradFBE;
                end
            end


        case 3 %truncated Newton method with conugate gradient
            % inactive (beta) and active (alpha) variables
            beta = find(abs(gs) <= lambda*gamma);
            alpha = setdiff(1:n, beta);
            result.dim(1,it) = nnz(alpha);
            % solve explicitly for inactive components
            d(beta) = -x(beta);
            % solve for active components
            cgForce = min(5e-1,(result.FPR(1,it)))*result.FPR(1,it);% superlinear
            cgReg = 1e-5*result.FPR(1,it);
            %             rhs = Atb(alpha)-lambda*sign(gs(alpha));
            %             [z,cgIter,~,how,exitflag] = conjGrad4TNreg(A(:,alpha),rhs,cgForce,cgReg,min(opt.cgMaxIter, result.dim(1,it)));
            %             d(alpha,1) = z-x(alpha);
            rhs2 = -(((A(:,alpha)*x(alpha,1)-b)'*A(:,alpha))'+lambda*sign(gs(alpha)));
            [d(alpha,1),cgIter,~,~,exitflag] = conjGrad4TNreg(A(:,alpha), ...
                rhs2,cgForce,cgReg,min(opt.cgMaxIter, result.dim(1,it)));
            cgIters(1,it) = cgIter;
            cgHow(1,it) = exitflag;
    end
    
    
    % line search (2 matvecs)
    Ad = A*d;
    AAd = At*Ad;
    Adres = Ad'*res;
    Ad2 = Ad'*Ad;
    
    
    
    % compute directional derivative
    slope = FPR'*(d-gamma*AAd);
    switch opt.linesearch
        case 1%Armijo
            tau   = 2;
            dorun = true;
            while dorun && tau > 1e-12
                tau = 0.5*tau;
                xd = x + tau*d;
                lossd  = loss + tau*Adres + 0.5*(tau^2)*Ad2;
                gradfd = gradf + tau*AAd;
                gsd = xd - gamma*gradfd;
                FBS = sign(gsd).*max(0,abs(gsd)-lambda*gamma);
                FPR = (xd-FBS)/gamma;
                % forward backward envelope in x+ss*d
                FBEd    = lossd + lambda*norm(FBS,1) - gamma*(gradfd'*FPR) + 0.5*gamma*(FPR'*FPR);
                dorun  = (FBEd > FBE+tau*1e-4*slope);
            end
        case 2%Wolfe
            [xd, tau, lossd, gradfd, gsd, FBEd,FBS, FPR] =...
                LemarechalLS(x, d, loss, Adres, Ad2, AAd, gradf, lambda, gamma, slope, FBE);
    end
    result.tau(1,it) = tau;
    
    
    % update iterate
    xp = x;
    x = xd;
    res = res+tau*Ad;
    loss = lossd;
    gradf = gradfd;
    gs = gsd;
    FBE = FBEd;
    
    
    % check for continuation
    if opt.cont && (lambda > lambda0) && norm(FPR) <= lambda*opt.tolc
        flagChangedLambda = true;
        icont = icont+1;
        if opt.direction == 1
            ls = sort(abs(gs)/gamma,1,'descend');
            lambda = max(lambda0, ls(m+1));
        else
            lambda = max(lambda0, opt.eta*lambda);
        end
        result.lambda(icont) = lambda;
        if icont==1
            result.its(icont) = it;
        else
            result.its(icont) = it-sum(result.its(1:icont-1))+1;
        end
        FBS = sign(gs).*max(abs(gs)-lambda*gamma,0);
        FPR = (x-FBS)/gamma;
        FBE = loss + lambda*norm(FBS,1) - gamma*(gradf'*FPR) + 0.5*gamma*norm(FPR,2)^2;
    end
    % Quick fix
    if opt.direction == 2
        gradFBEp = gradFBE;
        gradFBE = FPR/gamma-At*(A*FPR);
        gradFBEnrm = sqrt(gradFBE'*gradFBE);
    end
    
    
    alpha_prev = alpha;
    
end

%% ******************************************* PACK RESULTS
result.x     = FBS;
result.iterations = it-1;
result.pobj  = result.pobj(1:it);
result.times = result.times(1:it);
result.tau   = result.tau(1:it-1);
result.dim   = result.dim(1:it-1);
if opt.direction == 3,
    result.cgIters = cgIters(1:it-1);
end
result.res   = res;
result.loss  = loss;
result.gradf = gradf;
result.gamma = gamma;
result.gs    = gs;
result.L     = L;
result.alpha = alpha;
result.beta  = beta;
