function [p,k,res,how,exitflag] = conjGrad4TNreg(A,b,optTol,sigma,maxIter)
p = zeros(size(b));
r = b;
d = r;
rr = r'*r;
k = 0;
res = norm(r);
how = 'iterLimExc';
exitflag = 0;
while k < maxIter 
    Ad = A*d;
    AtAd = (Ad'*A)'+sigma*d;% form (A'*A+sigma*I)*d
    dAd = d'*AtAd;    
    % Conjugate Gradient
    alpha = rr/(dAd);
    p = p + alpha*d;
    r = r - alpha*AtAd;
    rr_new = r'*r;
    res = sqrt(rr_new);
    if  res <= optTol
        how = 'optTol';
        exitflag = 1;
        break
    end
    beta = rr_new/rr;
    d = r + beta*d;
    k = k + 1;    
    rr = rr_new;
end
