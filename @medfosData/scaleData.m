function scaleData( obj )
% SCALEDAT scales the data in the internal storage of the 
% the class and sets the flagf DATA_SCALED to TRUE
% Data should not be scaled more than once in order to avoid possible
% erranoeous behavior
% Scales only the data in training ang test sets! 

if obj.DATA_SCALED == 1
    display('Data is already scaled! Exiting...');
    return
end

for i = 1:obj.Ny
    vals = obj.Y(:,i);
    scaleFactor = max(abs(vals));
    
    obj.YTest(:,i) = obj.YTest(:,i) ./ scaleFactor;
    obj.YTrain(:,i) = obj.YTrain(:,i) ./ scaleFactor;
    
    YscaleMap(i) = scaleFactor;
end

for i = 1:obj.Nu
    vals = obj.U(:,i);
    scaleFactor = max(abs(vals));
    
    obj.UTest(:,i) = obj.UTest(:,i) ./ scaleFactor;
    obj.UTrain(:,i) = obj.UTrain(:,i) ./ scaleFactor;
    
    UscaleMap(i) = scaleFactor;
end


obj.DATA_SCALED = 1;
obj.YscaleMap = YscaleMap;
obj.UscaleMap = UscaleMap;

end

