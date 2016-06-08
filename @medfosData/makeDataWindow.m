function [ data , t, ynext, unext] = makeDataWindow(obj, stringName, currentIndex, currOutput)
%MAKEDATAWINDOW Reurns a structure of past values specified by currentIndex
% as specified by ARX polynomials structure

% There may be some issues when dealing with past Y values of different
% size

%currentIndex
inputInd = obj.uCross( currOutput, : );
outputInd = obj.yCross( currOutput, :);

Uwindow = zeros( max(inputInd), sum( inputInd > 0) );
Ywindow = zeros( max(outputInd), sum( outputInd > 0) );

if strcmp(stringName, 'Training')
    Y = obj.YTrain(:,outputInd > 0);
    U = obj.UTrain(:,inputInd > 0);
end

if strcmp(stringName, 'Test')
    Y = obj.YTest(:,outputInd > 0);
    U = obj.UTest(:,inputInd > 0);
end

count = 0;
for i=1:length(inputInd)
    
    if( inputInd(i) == 0 ), continue; end
    
    count = count + 1;
    nRegressors = inputInd( i );
    indexes = currentIndex: -1 : currentIndex - nRegressors + 1;
    
    Uwindow(1:length(indexes),count) = U(indexes, count);    
end

count = 0;
for i=1:length(outputInd)
    
    if( outputInd(i) == 0 ), continue; end
    
    count = count + 1;
    nRegressors = outputInd( i );
    indexes = currentIndex: -1 : currentIndex - nRegressors + 1;
    
    Ywindow(1:length(indexes),count) = Y(indexes, count);    
end

data.Y = Ywindow;
data.U = Uwindow;

if strcmp(stringName, 'Test')
    t = obj.YTest(currentIndex + 1 , currOutput );
    ynext = obj.YTest(currentIndex + 1, :);
    unext = obj.UTest(currentIndex + 1, :);
end

if strcmp(stringName, 'Training')
    t = obj.YTrain(currentIndex + 1 , currOutput );
    ynext = obj.YTrain(currentIndex + 1, :);
    unext = obj.UTrain(currentIndex + 1, :);
end
end

