function [ predict, predResiduals ] = createPredictModel( obj , methodStr, M )
%CREATEPREDICTMODEL creates predictor by optimizing a given criteria
%   methodStr == 'L1' uses sum of absolute values
%   methodStr == 'L2' uses squared 2 norm
%   methodStr == 'L1S' uses sum of absolute values with penalty on delta_Y
%   methodStr == 'L2S' uses squared 2 norm wint penalty on delta_Y
%   methodStr == 'huber' uses huber function
%   M is penalty parameter for delta_Y or huber function parameter
%TODO: ADD AN OPTION TO SOLVE LEAST SQUARES DIRECTLY 

%MOSEK seems to be the most robust one
ops = sdpsettings('solver','mosek','verbose',1); 

if strcmp(methodStr , 'L1')
    methodN = 1;
elseif strcmp(methodStr ,'L2')
    methodN = 2;
elseif strcmp(methodStr, 'huber') 
    methodN = 3;
elseif strcmp(methodStr ,'L2S')
    methodN = 4;
elseif strcmp(methodStr ,'L1S')
    methodN = 5;
else %default
    display('No method specifed, using L2 norm...');
    methodN = 2;
end
    
for i=1:obj.Ny
    RegressorsMatrix = obj.TrainMats{i}.Mat;
    TargetVector  = obj.TrainMats{i}.Targets;
    
    Constraints = [];
    Objective = 0;
    
    xStar = sdpvar( size(RegressorsMatrix, 2) , 1 );
    
    y_param = RegressorsMatrix*xStar;
    
    %L1 penalty
    if methodN == 1
        Objective = Objective + 1.0*norm( TargetVector - RegressorsMatrix*xStar , 1);
    %L2 penalty (least squares)
    elseif methodN == 2
        Objective = Objective + 0.5*norm( TargetVector - RegressorsMatrix*xStar , 2 )^2;
    %Huber penalty function
    elseif methodN == 3
        Objective = Objective + huber(TargetVector - RegressorsMatrix*xStar , M);
    %L2 penalty wit additonal penalty on Delat_Y weighted by penalty parameter 
    elseif methodN == 4
        Objective = Objective + 0.5*norm( TargetVector - y_param , 2 )^2 + M*norm( y_param(2:end) - y_param(1:end-1), 2 )^2;
    %L1 penalty wit additonal penalty on Delat_Y weighted by penalty parameter 
    elseif methodN == 5
        Objective = Objective + 0.5*norm( TargetVector - y_param , 1 ) + M*norm( y_param(2:end) - y_param(1:end-1), 1);    
    end
    
    sol = optimize( Constraints, Objective, ops );
    predict{i} = value(xStar);  

    %Not sure if this ever used. Remove?
    predResiduals{i} = TargetVector - RegressorsMatrix*predict{i};
end

end

