function [ Analytics, predictedTrajectory, measuredTrajectory] = predictAhead(  obj, predictorModel, K, lassoMap )
%PREDICTAHEAD Does a K steps ahead prediction on Test data. This means that
%Y values have been predicted for K steps back. U values have been used from the measured values  
%   Detailed explanation goes here

count = 0;

if nargin < 4 || isempty(lassoMap)
    for i=1:obj.Ny
        lassoMap{i} = logical( ones( length(predictorModel{i}), 1) );
    end
end


if isempty(obj.YscaleMap)
    scaleFactor = ones(obj.Ny,1);
else
    scaleFactor = obj.YscaleMap;
end


for k = obj.maxHist:length(obj.YTest) - 1 - K; 
    count = count + 1;
    
    for i = 1:obj.Ny
        [W,t,yNext, uNext] = makeDataWindow( obj,'Test', k, i);
        [h] = winToVec(obj, W, i);
        h = h(lassoMap{i});
        yPred(i) = h'*predictorModel{i}; %one step ahead prediction
        
        %save data only for update
        tempWin{i}.W = W;
        tempWin{i}.h = h;
        tempWin{i}.ypred = yPred;
    end
    
    Analytics{count}{1}.predictedY = yPred .* scaleFactor;
    Analytics{count}{1}.measuredY = yNext .* scaleFactor;
    
    % Predicted values
    tempWin2 = tempWin;
    
    %K steps ahead prediction loop
    for stepIterator = 1:K-1
        % iteratre over all given outputs
        for outIter=1:obj.Ny          
            tempWin{outIter}.W = updateDataWindow(obj, tempWin{outIter}.W, yPred, uNext, outIter );        
            h = winToVec(obj, tempWin{outIter}.W, outIter);
            h = h(lassoMap{outIter});
            yPTemp(outIter) = h'*predictorModel{outIter};
        end
        %get yNext and uNext that will be used to update current window
        [~,~, yNext, uNext] = makeDataWindow(obj,'Test',k+stepIterator,outIter);  
        yPred = yPTemp;
    
    %save predicted and measured values at each prediction length
    
    if ( size(yPred,1 ) ~= size( scaleFactor,1) )
        scaleFactor = scaleFactor';
    end
    
    %rescale the data
    Analytics{count}{1 + stepIterator}.predictedY = yPred .* scaleFactor;
    Analytics{count}{1 + stepIterator}.measuredY = yNext .* scaleFactor;
    
    end
    
    %save K steps ahead for easier plotting
    predictedTrajectory(count, :) = yPred .* scaleFactor;
    measuredTrajectory(count, :) = yNext .* scaleFactor ;    
    inputvalues(count,:) = uNext;

end

%figure, plot( predictedTrajectory ), hold on, plot( inputvalues(:,[1  4 7 ]) * 1000 );
%hold on, plot( measuredTrajectory ), hold on, plot( inputvalues(:,[2 5 8]) );
%hold on, plot( inputvalues(:,[10 11]) * 1100);

end


