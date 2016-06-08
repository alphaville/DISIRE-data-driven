function T = calcPerformanceIndices( obj, perfData , alpha)
%Calculates Average value at risk for given simulated data and coefficient
%alpha

Ny = obj.Ny;

predictionLen = length( perfData );
K = length( perfData{1} );

PMSE = zeros(predictionLen, Ny);
PRMSE = zeros(predictionLen, Ny);
NPRMSE = zeros(predictionLen,1);

for i=1:predictionLen
    summ = 0;
    
    for j=1:K
        errVec = perfData{i}{j}.predictedY - perfData{i}{j}.measuredY;
        summ = summ + dot(errVec, errVec); %2 squared norm of the error
    end
    
    PMSE(i,:) = 1/K * PMSE(i,:);
    PRMSE(i,:) = sqrt( PMSE(i,:) );
    
    %normed PRMSE beacause we have multiple outputs
    NPRMSE(i) = sqrt( 1 / K * summ ); 
end


T.NormedPRMSE = NPRMSE;

%calculate avar
figure;
h = histogram(NPRMSE,300);
values = h.Values;
            
pMass = sum(values);
p = values ./ pMass;
Z = 0:h.BinWidth:h.BinWidth*h.NumBins - h.BinWidth;

[a, mu] = obj.avar(Z, p, alpha);
T.AVAR = a;
return