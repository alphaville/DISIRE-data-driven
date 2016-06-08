function [ predict, predResiduals, lassoMap, sparInfo ] = createPredictModelLasso( obj, lambda )
%CREATEPREDICTMODEL Lasso regularization
%Solves a LASSO problem to find a subset of past input/output values
%resulting a "smaller" model
%Afterwards a new L2 penalty problem is solved to obtain predictor
%parameters
%MAYBE I CAN MAKE THIS MORE FLEXIBLE WITH RESPECT TO SECOND STAGE

ops = sdpsettings('solver','mosek','verbose',0); 
tol = 1e-8;

for i=1:obj.Ny
    RegressorsMatrix = obj.TrainMats{i}.Mat;
    TargetVector  = obj.TrainMats{i}.Targets;
    
    %YALMIP 
    Constraints = [];
    Objective1 = 0;
    Objective2 = 0;
    
    if (isempty(lambda) || (lambda == 0))
        lambda = 0.01 *  norm( RegressorsMatrix'*TargetVector, inf);
    end
    
    p = sdpvar( size(RegressorsMatrix, 2) , 1 );
    res = TargetVector - RegressorsMatrix*p;
    
    %LASSO 
    Objective1 = Objective1 + 0.5*(res'*res) +  lambda*norm(p,1);
    optimize(Constraints, Objective1,ops);
    
    %Determine a lasso map
    lassoMap{i} = abs(value(p)) > tol;
    spar = 1 - sum( lassoMap{i} ) / length( lassoMap{i} );
    sparInfo{i} = spar;
    
    %Second stage (debiasing)
    xStar = sdpvar( sum(lassoMap{i}), 1 );
    Objective2 = Objective2 + 0.5*norm( TargetVector - RegressorsMatrix(:,lassoMap{i}')*xStar,2)^2;
    optimize(Constraints, Objective2, ops);
    predict{i} = value( xStar );
   
    %not sure if this is used, but good for debugging
    predResiduals{i}  = TargetVector - RegressorsMatrix(:,lassoMap{i}')*predict{i};
end

end

