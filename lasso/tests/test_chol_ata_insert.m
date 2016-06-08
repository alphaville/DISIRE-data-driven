% GIST:
% https://gist.github.com/alphaville/4459f416c3d790b43502

! rm *~ ./tests/*~ 2> /dev/null

clear
repeat = 30;
N=120;
for r = 1:repeat,    
    for n=10:15,
        A = randn(N,n);
        L = chol(A'*A,'lower');
        c=rand(size(A,1), 1);
        for idx = 1:n+1,
            Ac = [A(:,1:idx-1) c A(:, idx:end) ];

            [l1, l2, p_tr, flag]=chol_ata_insert_col(L, A, c, idx);
            [~, p_] = sort(p_tr);
            assert(flag==0);

            % make sure the permutation is correct
            assert( norm([A c]-Ac(:,p_)) < 1e-14  );

            I=eye(n+1);
            Pmat = I(:, p_tr);
            L_updated = [L zeros(size(L,1),1); l1' l2];
            T = L_updated*L_updated';
            assert(  norm(T-[A c]'*[A c]) < 1e-10  );
            assert(  norm(Ac'*Ac - Pmat'*T*Pmat) < 1e-10  );
            assert(  norm(Ac'*Ac - T(p_tr, p_tr)) < 1e-10  );
        end
    end
end