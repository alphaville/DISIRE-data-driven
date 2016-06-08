% solve a lasso problem using ForBES

%close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 280; n = 5000;
x_orig = sprandn(n, 1, 20/n); % generate random sparse model
A = sprandn(m, n, 50/n); % generate random sparse design matrix
b = A*x_orig + randn(m, 1)/10; % compute labels and add noise

fprintf('%d nonzero features\n', nnz(A));
fprintf('%.2f nnz per row\n', nnz(A)/numel(A)*n);

% for lam >= lam_max the solution is zero

lam_max = norm(A'*b,'inf');


Hess = @(x) A'*(A*x);eigsOpt.issym = 1;eigsOpt.tol = 1e-5;
Lf = eigs(Hess, n, 1, 'LM', eigsOpt);
gamma = 0.999/Lf;

% run some methods

%

lam = 0.6*lam_max;


clear opt2
fprintf('\nFBN-CG\n');
opt2.maxit = 100;
opt2.direction = 3
opt2.linesearch = 1;
opt2.cont = false;

opt2.tol=1e-8;
opt2.tolc = 1e-4;
opt2.gamma=gamma;
opt2.chol_update = false;



tic;out_FBN = FBNewtonL1LS(A, b, lam,opt2);toc
fprintf('\n');
fprintf('iterations   : %d\n', out_FBN.iterations);
fprintf('time         : %7.4e\n', out_FBN.times(end));
fprintf('residual     : %7.4e\n', out_FBN.FPR(end));
plot(out_FBN.times, out_FBN.pobj,'-o');
grid on
%%
fprintf('\n L-BFGS\n');
opt2.direction = 2;
opt2.linesearch = 2;
tic;out_LBFGS = FBNewtonL1LS(A, b, lam,opt2);toc
fprintf('\n');
fprintf('iterations   : %d\n', out_LBFGS.iterations);
fprintf('time         : %7.4e\n', out_LBFGS.times(end));
fprintf('residual     : %7.4e\n', out_LBFGS.FPR(end));
fprintf('diff from #1 : %7.4e\n', norm(out_FBN.x-out_LBFGS.x, Inf));


%%
f = quadLoss(1, zeros(m,1));
aff = {A, -b};
g = l1Norm(lam);
x0 = zeros(n, 1);

opt.maxit = 10000;
opt.tol = 1e-10;
opt.adaptive = 0;
opt.display = 1;
opt.memory = 5;
% run some methods
fprintf('\nFBS\n');
opt_fbs = opt;
opt_fbs.method = 'fbs';
opt_fbs.variant = 'fast';
tic; out_fbs = forbes(f, g, x0, aff, [], opt_fbs); toc
fprintf('\n');
fprintf('message    : %s\n', out_fbs.message);
fprintf('iterations : %d\n', out_fbs.iterations);
fprintf('matvecs    : %d\n', out_fbs.operations.C1);
fprintf('prox       : %d\n', out_fbs.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.residual(end));

%%
fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
opt_lbfgs.variant = 'basic';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('\n');
fprintf('message    : %s\n', out_lbfgs.message);
fprintf('iterations : %d\n', out_lbfgs.iterations);
fprintf('matvecs    : %d\n', out_lbfgs.operations.C1);
fprintf('prox       : %d\n', out_lbfgs.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.residual(end));

assert( norm(out_FBN.x-out_lbfgs.x, Inf) < 1e-7 );



%%

randn('seed', 0);
rand('seed',0);

m = 1500;       % number of examples
n = 5000;       % number of features
p = 100/n;      % sparsity density

x0 = sprandn(n,1,p);
A = randn(m,n);
A = A*spdiags(1./sqrt(sum(A.^2))',0, n, n); % normalize columns
b = A*x0 + sqrt(0.001)*randn(m,1);

lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max;

[x, history] = lasso_ADMM_boyd(A, b, lambda, 1.0, 1.0);


%% Test QP
Q = [ 3       1       1       0
    1       3       0      -1
    1       0       3       0
    0      -1       0       5];
eig(Q)
assert(all(eig(Q)>0))
q=[1;2;3;4];
A=[ 1       1       0       0
    0       0       1       2];
b=[1;0];
H = [1 1 1 1];

x = sdpvar(4,1);
J = 1/2*x'*Q*x + q'*x;
solvesdp([-1 <= H*x <= 1, A*x == b],J)

x_star = double(x)


%% PI approximation
X=[];
tol = 1e-8;
e=exp(1);
for p=4:50000,
    % p/q
    for q=ceil(p/3):floor(p/2),
        if (abs(p-q*e)<q*tol),
            fprintf('%g/%g\n', p, q);
        end
    end
end


%% Contour plot of myfun0
x=[]; y=[]; z=[];
for i=-1:0.01:1,
    for j=-1:0.01:1,
        x=[x; i];
        y=[y; j];
        z=[z;myfun01(i,j)];        
    end
end