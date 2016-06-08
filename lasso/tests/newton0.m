%% TEST beta_0, beta_1, beta_2
clc;
clear all;

% N = 50; n = 1500;
% 
% A = randn(N, n);
% y = A*ones(n,1);
% lambda=3000;
% x0 = 0.1*rand(n,1);

m = 20;
n = 2000;
mstar = 10;
lambda = 5;
[A,b,x,y] = generateL1(m,n,mstar,lambda);
x0 = zeros(n,1);
A = full(A);

%%
% Solve with Newton
clc
A = randn(10,100);
H = diag(.9./diag(A'*A));
A = A*H;
y = 10*rand(10,1);
ops.max_iter = 600;
ops.tol = 1e-7;
ops.sigma = 1e-4;
x0 = rand(100,1);
lambda=0.1;
[x_star, details ] = lasso_newton(A, y, lambda, x0, ops);
subplot(211)
semilogy(details.cost,'-o')
subplot(212)
semilogy(details.tau,'-o')




%% Sovle with YALMIP
x = sdpvar(size(x0,1),1);
J = 0.5*(A*x-y)'*(A*x-y) + lambda*norm(x,1);
ops = sdpsettings;
ops.solver='mosek'; % try also with sedumi
ops.verbose=0;
o=solvesdp([], J, ops); 
x = double(x);
J = double(J);
err = (0.5*(A*x_star-y)'*(A*x_star-y) + lambda*norm(x_star,1) - J)/J;
fprintf('Error : %2.3f %% \n', err*100);
close; figure;
stem(x,'-o');
hold on;
stem(x_star,'Linewidth',2); legend('YALMIP','NEWTON');
axis tight;
o


%%
E = inf(8,8);
V = inf(8,8);
[E,V]=assignme(E,V,1,2,5,10);
[E,V]=assignme(E,V,2,3,6,8);
[E,V]=assignme(E,V,2,6,4,1);
[E,V]=assignme(E,V,6,7,9,1);
[E,V]=assignme(E,V,6,3,2,1);
[E,V]=assignme(E,V,4,7,3,2);
[E,V]=assignme(E,V,3,4,1,5);
[E,V]=assignme(E,V,3,6,1,0.5);
[E,V]=assignme(E,V,4,5,10,10);
[E,V]=assignme(E,V,1,8,4.5,0.5);
[E,V]=assignme(E,V,7,5,7,2);
[E,V]=assignme(E,V,8,5,23,2);
[E,V]=assignme(E,V,8,3,9,2);
% V: degrees of freedom of chi squared
clc


path = [1 8 3 4 5];
Mean=0;
Var = 0;
lambda = 10;
for i=1:(length(path)-1),
    Mean = Mean + E(path(i),path(i+1));
    Var = Var + V(path(i),path(i+1));
end
fprintf('Mean = %g, Var = %g, Cost = %g\n', Mean, Var, (Mean+lambda*Var))
fprintf('P[X > 30] = %g%%\n', 100*(1 - normcdf(30,Mean, sqrt(Var))))



%% 

% solve a lasso problem using ForBES

close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 100; % number of observations
n = 500; % number of features
x_orig = sprandn(n, 1, 20/n); % generate random sparse model
A = sprandn(m, n, 50/n); % generate random sparse design matrix
b = A*x_orig + randn(m, 1)/10; % compute labels and add noise

fprintf('%d nonzero features\n', nnz(A));
fprintf('%.2f nnz per row\n', nnz(A)/numel(A)*n);

% for lam >= lam_max the solution is zero

lam_max = norm(A'*b,'inf');
lam = 0.0*lam_max;

f = quadLoss(1, zeros(m,1));
aff = {A, -b};
g = l1Norm(lam);
x0 = zeros(n, 1);

% set common options

% run some methods
fprintf('\nFBN\n');
opt.cont = false;
tic;out_FBN = FBNewtonL1LS(A, b, lam, opt);toc
fprintf('\n');
fprintf('iterations : %d\n', out_FBN.iterations);
fprintf('time       : %7.4e\n', out_FBN.times(end));
fprintf('residual   : %7.4e\n', out_FBN.FPR(end));


%%


opt.maxit = 10000;
opt.tol = 1e-10;
opt.adaptive = 0;
opt.display = 1;
opt.memory = 5;

fprintf('\nFBS\n');
opt_fbs = opt;
opt_fbs.method = 'fbs';
opt_fbs.variant = 'fast';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('\n');
fprintf('message    : %s\n', out_fbs.message);
fprintf('iterations : %d\n', out_fbs.iterations);
fprintf('matvecs    : %d\n', out_fbs.operations.C1);
fprintf('prox       : %d\n', out_fbs.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
opt_lbfgs.variant = 'global';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('\n');
fprintf('message    : %s\n', out_lbfgs.message);
fprintf('iterations : %d\n', out_lbfgs.iterations);
fprintf('matvecs    : %d\n', out_lbfgs.operations.C1);
fprintf('prox       : %d\n', out_lbfgs.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.residual(end));


norm(out_FBN.x - out_fbs.x, Inf)

