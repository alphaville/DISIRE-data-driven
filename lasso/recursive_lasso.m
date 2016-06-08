function [times, iters] = recursive_lasso(n, f_sparse, sigma, options)
%
% n             window size
% f_sparse      density of the data stream [0~1]
% sigma         variance of the observation noise (warning: this is often
%               denoted as sigma^2)
% options       optional structure with options;

% take care of optional options
if nargin < 4 || ~isfield(options, 'fbs') || isempty(options.fbs),options.fbs=false; end
if nargin < 4 || ~isfield(options, 'lbfgs') || isempty(options.lbfgs),options.lbfgs=false; end
if nargin < 4 || ~isfield(options, 'admm') || isempty(options.admm),options.admm=false; end
if nargin < 4 || ~isfield(options, 'l1ls') || isempty(options.l1ls),options.l1ls=false; end
if nargin < 4 || ~isfield(options, 'Nsim') || isempty(options.Nsim),options.Nsim=20; end


s = RandStream('mt19937ar', 'Seed', 'shuffle');
RandStream.setGlobalStream(s);
reset(s);

%f_sparse = 0.1; % sparsity

%n = 5000; % number of observations
m = 4*ceil(f_sparse*n); % number of features


assert(m<n)

x_orig = sprandn(n, 1, f_sparse); % generate random sparse model
A = randn(m, n)/sqrt(m); % generate random sparse design matrix

b = A*x_orig + randn(m, 1)/10; % compute labels and add noise
% lam_max = norm(A'*b,'inf');
% lam = 0.25*lam_max;

Hess = @(x) A'*(A*x);
eigsOpt.issym = 1;
eigsOpt.tol = 1e-5;
Lf = eigs(Hess, n, 1, 'LM', eigsOpt);
gamma = 0.999/Lf;
options.maxit = 300;
options.cont = 1;
options.tol=1e-8;
options.tolc = 1e-8;
options.gamma=gamma;

%% Run algorithms with streaming data

%sigma = 0.005; % noise

% Construct long stream of data
N = 1e6;                                          % stream length
N_sim = min(options.Nsim, floor(N/n));            % simulation length
X=zeros(N,1); idx = randperm(N); nnz_X = floor(f_sparse*N);
idx_nonzero = idx(1:nnz_X);
idx_minus = idx_nonzero(binornd(1, 0.5, [nnz_X,1])==1);
idx_plus  = setdiff(idx(1:nnz_X), idx_minus);     % indices 
magic_val= 8 * sqrt(2*sigma*log(N_sim));
X(idx_plus) =  ( 1+rand(length(idx_plus), 1))*magic_val + 1e-9;   
X(idx_minus) = (-1-rand(length(idx_minus),1))*magic_val - 1e-9; 


A_i = A;
Axi = [];               % we store and update the values of A*x_i to accelerate the algorithm
options.x = zeros(n,1); % initial guess (zero)


iters = zeros(N_sim,2); % number of iterations
times = zeros(N_sim,5); % times = [FBN, FBS, LBFGS, ADMM, Gurobi]

for i=1:N_sim,
    % Sample from the stream - use window of length n
    x_i = X(i:i+n-1);
    
    % Observe output
    y_i = A_i * x_i + sigma*randn(m,1);
    
    % Set lambda
    lam = 4  * sqrt(2 * sigma * log(n));
    
    % [FBN - Our algorithm]
    if ~isempty(Axi), out.res = Axi - y_i; end
    if i==1,
        options.direction = 3; 
    else
        options.direction = 1; 
    end
    options.linesearch = 1; 
    opt.chol_update=true;
    opt.cont = true;
    start=tic;
        out = FBNewtonL1LS(A_i, y_i, lam, options);
    time_fbn=toc(start);
    iters(i,1) = out.iterations; 
    times(i,1) = time_fbn;
    fprintf('%d\t FBN iterations : %g, \t FBN time: %f \n', i, out.iterations, time_fbn);
    Axi = out.res + y_i; % keep track of A_i*x_i to speed things up
    
    if i>1,
        
        % [FBS with forbes]
        if options.fbs,
            f = quadLoss(1, zeros(m,1)); aff = {A_i, -y_i}; g = l1Norm(lam);
            x0 = options.x;
            opt_fbs.maxit = 10000;    opt_fbs.tol = options.tol;
            opt_fbs.adaptive = 0;     opt_fbs.display = 0;
            opt_fbs.method = 'fbs';   opt_fbs.variant = 'fast';
            start_fbs = tic;
            out_fbs = forbes(f, g, x0, aff, [], opt_fbs);
            times(i, 2) = toc(start_fbs);
            fprintf('%d\t FBS iterations : %g, \t FBS time: %f \n', i, out_fbs.iterations, times(i, 2));
        end
        
        % [ADMM Boyd]
        if options.admm,
            start_admm=tic; x_admm = lasso_ADMM_boyd(A_i, y_i, lam, 1.0, 1.0);
            times(i, 3) = toc(start_admm);
            fprintf('%d\t ADMM time : %f \n', i, times(i, 3));
        end
        
        % [LBFGS with forbes]
        if options.lbfgs,
            f = quadLoss(1, zeros(m,1)); aff = {A_i, -y_i}; g = l1Norm(lam);
            x0 = options.x;
            opt_lbfgs.maxit = 10000;    
            opt_lbfgs.tol = options.tol;
            opt_lbfgs.adaptive = 0;     
            opt_lbfgs.display = 0;          
            opt_lbfgs.memory = 5;
            opt_lbfgs.method = 'lbfgs';
            opt_lbfgs.variant = 'basic';
            start_lbfgs = tic;
            out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
            times(i, 4) = toc(start_lbfgs);
            fprintf('%d\t LBFGS iterations : %g,  LBFGS time: %f \n', i, out_lbfgs.iterations, times(i, 3));
        end
        
        
        % [L1LS] -- see https://stanford.edu/~boyd/l1_ls/
        if options.l1ls,
            start_l1_ls = tic; [x,status] = l1_ls(A_i, y_i, lam); times(i, 5) = toc(start_l1_ls);
            fprintf('%d\t L1LS time: %f \n', i, times(i, 5));
        end
    end
    
    % plot
    stairs(circshift(out.x,-1)); axis([0 n -3 3]);hold off;pause(0.01);
    
    % Permute A_i
    A_i = circshift(A_i', -1)';
    
    % Update the initial guess
    options.x = circshift(out.x,-1);
    
    
end

%bar(mean(times(2:end,:)))
%Labels = {'FBN', 'FBS', 'LBFGS', 'ADMM', 'L1LS'};
%set(gca, 'XTick', 1:5, 'XTickLabel', Labels);



