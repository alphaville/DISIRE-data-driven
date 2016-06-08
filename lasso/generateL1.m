function [A,b,x,y] = generateL1(m,n,mstar,rho)
A = sparse(m,n);
B = sprand(m,n,5e-2);
[i,j]=find(B);
B = sparse(i,j,-1,m,n)+2*B;
v = rand(m,1);
y = v/norm(v);
[Bys,is] = sort(abs(B'*y),'descend');
B = B(:,is);
for i=1:n
    if i<=mstar
        alpha=1/(B(:,i)'*y);
    elseif abs(B(:,i)'*y)<=0.1 && i>mstar
        alpha = 1;
    else
        alpha = rand/abs(B(:,i)'*y);
    end
    A(:,i)=alpha*B(:,i);
end
x=zeros(n,1);
for i=1:mstar
        xi = rho/sqrt(mstar)*rand;
        x(i,1) = xi*sign(A(:,i)'*y);
end
b = y+A*x;