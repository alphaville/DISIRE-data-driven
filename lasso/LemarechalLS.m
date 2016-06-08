function [xd, tau, lossd, gradfd, gsd, FBEd,FBS, FPR, exitflag] = LemarechalLS(x, d, loss, Adres, Ad2, AAd, gradf, lambda, gamma, slope, FBE)
lsopt.sigma = 0.9;
lsopt.nbracket = 100;
lsopt.interp = 1;
lsopt.delta = 0.1;
lsopt.rho = 5;
lsopt.theta = 0.49;
lsopt.progTol = 0;

tau = 1;
wolfe_hi = lsopt.delta*slope;
wolfe_lo = lsopt.sigma*slope;
a = 0; fa = FBE; dfa = slope;
tprev = a; fprev = fa; dfprev = dfa;
b = inf; % upper bound
rho = lsopt.rho;
theta = lsopt.theta;
exitflag = 1;
for it = 1:lsopt.nbracket
    xd = x + tau*d;
    lossd  = loss + tau*Adres + 0.5*(tau^2)*Ad2;
    gradfd = gradf + tau*AAd;
    gsd = xd - gamma*gradfd;
    FBS = sign(gsd).*max(0,abs(gsd)-lambda*gamma);
    FPR = (xd-FBS)/gamma;
    % forward backward envelope in x+ss*d
    FBEd    = lossd + lambda*norm(FBS,1) - gamma*(gradfd'*FPR) + 0.5*gamma*(FPR'*FPR);
    if FBEd > FBE + tau*wolfe_hi
        b = tau; fb = FBEd;
        if lsopt.interp == 1
            tn = LemarechalQuadInterp(a,fa,dfa,b,fb);
            % safeguard
            tn = min(tn,b - theta*(b - a));
            tn = max(tn,a + theta*(b - a));
        else
            tn = 0.5*(a + b);
        end
        tau = tn;
    else
        dFBE = FPR'*(d-gamma*AAd);
        if dFBE < wolfe_lo
            a = tau; fa = FBEd; dfa = dFBE;
            if b == inf
                % extrapolate
                if lsopt.interp% we always have dfprev
                    tn = LemarechalCubInterp(tprev,fprev,dfprev,a,fa,dfa);
                    % safeguard
                    tn = max(tn,rho*tprev);
                else
                    tn = rho*tau;
                end
            else
                % interpolate
                if lsopt.interp == 1
                    tn = LemarechalQuadInterp(a,fa,dfa,b,fb);
                    % safeguard
                    tn = min(tn,b - theta*(b - a));
                    tn = max(tn,a + theta*(b - a));
                else
                    tn = 0.5*(a + b);
                end
            end
            tprev = tau;fprev = fa;dfprev = dfa;
            tau = tn;
        else
            exitflag = 0;
            break;
        end
    end
    
    if (b-a) <= lsopt.progTol
        exitflag = 2;
        break;
    end
end

end

function t = LemarechalQuadInterp(t0,f0,df0,t1,f1)
% Minimizer of interpolant belongs to [0,t1]
q = f1-f0-t1*df0;
q = 2*(f1-f0-(t1-t0)*df0)/(t1-t0)^2;
if q > 0%quadratic is strongly convex
    c2 = df0-t0*q;
    t = -c2/q;
else
    t = -1;
end
end

function t = LemarechalCubInterp(t1,f1,df1,t2,f2,df2)
% t1, t2 might not be sorted
delta = t2 - t1 ;
if  delta == 0
    t = t1;
else
    d1 = df1 + df2 - 3*(f2-f1)/delta;
    d2 = d1^2 - df1*df2 ;
    if  d2 < 0 % /* complex roots, use secant method */
        if ( abs(df1) < abs (df2) )
            t = t1 - (t1 - t2)*(df1/(df1-df2)) ;
        elseif ( df1 ~= df1 )
            t = t2 - (t1 - t2)*(df2/(df1-df2)) ;
        else
            t = - 1 ;
        end
    else
        % first way: from Hager-Zhang code
        %                 d2 = sqrt(d2)*sign(delta);
        %                 v1 = df1 + d1 - d2 ;
        %                 v2 = df2 + d1 + d2 ;
        %                 if ( (v1 == 0) && (v2 == 0) )
        %                     t = -1;
        %                 elseif ( abs (v1) >= abs (v2) )
        %                     t = t1 + delta*df1/v1 ;
        %                 else
        %                     t = t2 - delta*df2/v2 ;
        %                 end
        %
        % second way: from Bonnans, Lemarechal
        d2 = sqrt(d2)*sign(delta);
        v1 = df2 + d2 - d1;
        v2 = df2 - df1 + 2*d2;
        if ( (v1 == 0) && (v2 == 0) )
            t = -1;
        else
            t = t2 - delta*(v1/v2);
        end
    end
end
end

