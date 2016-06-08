classdef medfosData < handle
    %A placeholder for mefos data and various methods for system ident
    
    properties (GetAccess=public, SetAccess=private)
        %basic data about of the experiment
        YTrain = [];
        UTrain = [];
        YTest  = [];
        UTest  = [];
        Y = [];
        U = [];
        
        %number of the input and output poperties used
        Ny = 0;
        Nu = 0;
        
        % description of ARX model polynomials
        uCross = [];
        yCross = [];
        
        %Various counters
        maxHist = 0;
        testCounter = 0;
        
        %Training and Testing data in matrix form
        TrainMats = [];
        TestMats = []
        
        %Various flags
        DATA_SCALED = 0;
        
        %Various maps
        YscaleMap = [];
        UscaleMap = [];
        
    end
    
    methods
        function self = medfosData( dataStruct, targetNames, inputNames, dataPartition, polyStruct ) 
            tic
            
            Ny = length(targetNames);
            Nu = length(inputNames);

            k1 = dataPartition.Test.start;
            k2 = dataPartition.Test.end;

            inLen  = max(size(inputNames));
            outLen = max(size(targetNames));

            Y = [];
            U = [];

            %extract data from structure
            for i=1:inLen
                U = [U, dataStruct.(inputNames{i}) ];
            end

            for i=1:outLen
                Y = [Y, dataStruct.(targetNames{i}) ];
            end
            %extract Test data
            k1 = dataPartition.Training.start;
            k2 = dataPartition.Training.end;
            count = 0;
            
            for i = k1:k2
                count = count + 1;
                YTrain(count,:) = Y(i,:);
                UTrain(count,:) = U(i,:);
            end
            %extract Training data
            k1 = dataPartition.Test.start;
            k2 = dataPartition.Test.end;
            count = 0;
            
            for i = k1:k2
                count = count + 1;
                YTest(count,:) = Y(i,:);
                UTest(count,:) = U(i,:);
            end
            
            mt1 = max( max( polyStruct.inCross ) );
            mt2 = max( max( polyStruct.outCross) );
            
            maxHist = max( mt1, mt2 );
            
            %set properties
            self.Y = Y;
            self.U = U;
            self.YTrain = YTrain;
            self.UTrain = UTrain;
            self.YTest = YTest;
            self.UTest = UTest;
            self.Ny = Ny;
            self.Nu = Nu;
            self.uCross = polyStruct.inCross;
            self.yCross = polyStruct.outCross;
            self.maxHist = maxHist;
            toc
        end
        
        %other methods
        vec = winToVec( obj, W, indY );
        createTrainingMat(obj);
        scaleData(obj);
        [ predict, predResiduals ] = createPredictModel( obj, methodStr, M);
        [ predictor, residuals ] = createDummyPredictor(obj, N);
        [ data , t, ynext, unext] = makeDataWindow(obj, stringName, currentIndex, currOutput);
        newWin = updateDataWindow(obj,W, newY, newU, output);
        T = calcPerformanceIndices( obj, perfData , alpha);
        [a,mu] = avar( obj, Z, p, alpha)
        BFR = calculateBFR(obj,modelPrediction, dummyPrediction, measuredValues);
        poles = calculatePoles(obj, predictor )
    end
    
end
