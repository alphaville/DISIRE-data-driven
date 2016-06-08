! rm *~ ./tests/*~ 2> /dev/null
clear
A=rand(100,5);
AtA = A'*A;


N = 300;
n = 9;


% Add and remove
clear
N = 800;
n = 15;

A = rand(N,n);


alpha     = [    2  3  5  7  8  11    ];
alpha_bar = [1   2     5  7     11  12];
L = chol(A(:, alpha)'*A(:, alpha), 'lower');


[alpha_to_kill, ikill] = setdiff(alpha, alpha_bar);
alpha_copy = alpha;
alpha_copy(ikill) = [];


[L11bar, L31bar, L33bar] = chol_ata_remove_col(L, ikill);
L_kill = [ L11bar  zeros(size(L11bar,1), size(L33bar,2))
    L31bar  L33bar ];
A_kill = A(:,alpha_copy);


L = L_kill;
idx_to_add = setdiff(alpha_bar, alpha);
for i=idx_to_add,    
    [l1,l2, flag]=chol_ata_append_col(L, A_kill, A(:, i));
    assert(flag==0);
    L = [L zeros(size(L,1), 1)
        l1' l2];
    A_kill = [A_kill A(:,i)];
end


%% CHOL_ATA_UPDATE
clear
N = 250;
n = 15;
A = rand(N,n);
alpha     = int32(1:10);
alpha_bar = int32(1:11);
L = chol(A(:, alpha)'*A(:, alpha), 'lower');


[L_, perm_, out] = chol_ata_update(A, L, alpha, alpha_bar);

assert(length(perm_)==length(alpha_bar), 'Bug 27')

assert( norm(L_*L_' - A(:, perm_)'*A(:, perm_), Inf)  < 1e-10 , 'Wrong result');


%% EMPTY AND FILL IN AGAIN

N = 250;
n = 15;
A = rand(N,n);
alpha     = int32([    2  3  5  7  8  11    ]);
alpha_bar = int32([2  1   5  11  7 12  14]);
L = chol(A(:, alpha)'*A(:, alpha), 'lower');
[L_, perm_, out] = chol_ata_update(A, L, alpha, alpha_bar);
assert( norm(L_*L_' - A(:, perm_)'*A(:, perm_), Inf)  < 1e-10);


%% MORE COLUMNS
clear
N = 10;
n = 100;
A = rand(N,n);
alpha=1:10:100;
alpha_bar=1:2:100;


L = chol(A(:, alpha)'*A(:, alpha), 'lower');

[L_, perm_, flag] = chol_ata_update(A, L, alpha, alpha_bar);
