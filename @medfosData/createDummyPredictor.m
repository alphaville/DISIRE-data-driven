function [ predictor, residuals ] = createDummyPredictor( obj, N )
%Creates a dummy predictor that averages last N available values
%used to compare with more sophisticated predictors
%If N > any of the vales on obj.yCross values it will be replaced with
%aforementioned value

if N <= 0
    error('N has to be a positive integer') %TODO test for integer
end

for i=1:obj.Ny
    if( N > obj.yCross(i,i) )
        N = obj.yCross(i,i);
        display('N truncated to fit past output values history size');
    end
end

for i=1:obj.Ny
    RegressorsMatrix = obj.TrainMats{i}.Mat;
    TargetVector  = obj.TrainMats{i}.Targets;
    
    predictorSize = sum( obj.yCross(i,:) ) + sum( obj.uCross(i,:) );
    predictorVector = zeros( predictorSize, 1 );
    
    startIndex = sum( obj.yCross(i,1:i-1) ) + 1;
    endIndex = startIndex + N - 1;
    predictorVector(startIndex:endIndex) = 1/N;
    
    predictor{i} = predictorVector;
    residuals{i} = TargetVector - RegressorsMatrix*predictor{i};
end

end

