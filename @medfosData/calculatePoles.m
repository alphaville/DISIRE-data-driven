function poles = calculatePoles(obj, predictor )

for i=1:obj.Ny
    startIndex = sum( obj.yCross(i,1:i-1) ) + 1;
    endIndex = startIndex + obj.yCross(i,i) - 1;
    poles{i} = roots([1; -predictor{i}(startIndex:endIndex)]);
end

end
