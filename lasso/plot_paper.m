%% Experiment #1
% sparsity : f_sparse = 0.1
% experiment with n
set(0,'DefaultAxesFontSize',14)
n = [1000;
    5000;
    10000;
    20000];
times = [   0.0118    0.3603    0.0350    0.2657    0.1905
            0.1724    4.7645    1.5174    2.8500    2.2100
            0.6984   18.2312    9.5608   10.5319    8.2165
            3.3904   76.3005   63.0690   42.6574   31.6252];

figure(1);
times=times(:,[1 2 4 5]);
loglog(n, times(:,1),'-s', 'LineWidth',2, 'MarkerSize', 8); hold on
loglog(n, times(:,2),'-d', 'LineWidth',2, 'MarkerSize', 8)
loglog(n, times(:,3),'-o', 'LineWidth',2, 'MarkerSize', 8)
loglog(n, times(:,4),'-x', 'LineWidth',2, 'MarkerSize', 8)
axis tight
legend('FBN', 'FISTA', 'ADMM', 'L1LS','Location','SouthEast')
xlabel('Window size')
ylabel('Average runtime [s]')
grid on
hold off;

%%
clear
set(0,'DefaultAxesFontSize',14)
figure(2);
load runtimes_sparsity.mat
f_sparsity = f_sparsity *100;
results = results(:, [1 2 4 5])
semilogy(f_sparsity, results(:,1),'-s', 'LineWidth',2, 'MarkerSize', 8)
hold on
semilogy(f_sparsity, results(:,2),'-d', 'LineWidth',2, 'MarkerSize', 8)
semilogy(f_sparsity, results(:,3),'-o', 'LineWidth',2, 'MarkerSize', 8)
semilogy(f_sparsity, results(:,4),'-x', 'LineWidth',2, 'MarkerSize', 8)
axis([0 15 0 5])
legend('FBN', 'FISTA', 'ADMM', 'L1LS','Location','SouthEast')
xlabel('Sparsity [%]')
ylabel('Average runtime [s]')
grid on
hold off;


%% Experiment #3 (Sparsity)
clear ops
n=7000;
sigma=0.1;
f_sparsity=[0.05];

ops.Nsim = 50;

results = [];
for f=f_sparsity,
    fprintf('Running for f = %g\n', f);
    times = recursive_lasso(n, f, sigma, ops);
    times = mean(times(2:end,:))
    results = [results; times];
end


%% Experiment #4
clear
n = 10000;
f_sparse = 0.2;


ops.Nsim = 25;
ops.l1ls = 0;
ops.tol = 1e-6;


[times, iters]   = recursive_lasso(n, f_sparse, 0.1, ops)


%% Experiment #5 (Sigma)
clear ops
n=7000;
sigmas=[0.001 0.01 0.1 0.2 0.25 0.5 1.0 1.5];
f=0.05;

ops.Nsim = 40;
ops.fbs = true;
ops.lbfgs = true;
ops.admm = true;
ops.l1ls = true;

results = [];
for sigma=sigmas,
    fprintf('Running for sigma = %g\n', sigma);
    times = recursive_lasso(n, f, sigma, ops);
    times = mean(times(2:end,:))
    results = [results; times];
end


