function BFR = calculateBFR(obj,modelPrediction, dummyPrediction, measuredValues)
%Calculates how good our model is against a very basic one
%Prefered values are close to 1

rModel = modelPrediction - measuredValues;
rDummy = dummyPrediction - measuredValues;

if size(rModel,2) > size(rModel,1)
    rModel = rModel';
end

if size(rDummy,2) > size(rDummy,1)
    rDummy = rDummy';
end

for i=1:obj.Ny
    BFR(i) = 1 - (rModel(:,i)'*rModel(:,i)) / (rDummy(:,i)'*rDummy(:,i));
end

return
