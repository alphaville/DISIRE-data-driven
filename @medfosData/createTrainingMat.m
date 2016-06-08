function  createTrainingMat(obj)
%CREATETRAININGMAT creates a big matrix which is to be used for model
%training. It basically stacks row by row all the values used to predict
%current target value
%Output is cell of Ny matrices written directly to class itself


retMat = cell(obj.Ny, 1);
totalLen = length(obj.YTrain) - obj.maxHist;

%preallocate for speed
for i=1:obj.Ny
    totalWidth = sum(obj.yCross(i,:)) + sum(obj.uCross(i,:));
    retMat{i}.Mat = zeros( totalLen, totalWidth );
    retMat{i}.Targets = zeros( totalLen, 1);
end

for i=1:obj.Ny
    count = 0;
    
    for k=obj.maxHist:length(obj.YTrain) - 1
        [W,t,~,~] = makeDataWindow( obj,'Training', k, i );
        h = winToVec(obj, W, i); %matvec operation
        
        count = count + 1;
        retMat{i}.Mat(count,:) = h';
        retMat{i}.Targets(count,1) =  t;
    end

        
       
end

obj.TrainMats = retMat;

retMat = cell(obj.Ny, 1);
totalLen = length(obj.YTest) - obj.maxHist + 1;

for i=1:obj.Ny
    totalWidth = sum(obj.yCross(i,:)) + sum(obj.uCross(i,:));
    retMat{i}.Mat = zeros( totalLen, totalWidth );
    retMat{i}.Targets = zeros( totalLen, 1);
end

for i=1:obj.Ny
    count = 0;
    
    for k=obj.maxHist:length(obj.YTest) - 1
        [W,t,~,~] = makeDataWindow( obj,'Test', k, i );
        h = winToVec(obj, W, i); %matvec operation
        
        count = count + 1;
        retMat{i}.Mat(count,:) = h';
        retMat{i}.Targets(count,1) =  t;
    end

        
       
end
obj.TestMats = retMat;

end
