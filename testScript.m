clear all;
load data/medfos

targetNms{1}  = 'WBF_Z01_ZoneTempHSIMV';
targetNms{2}  = 'WBF_Z02_ZoneTempHSIMV';
targetNms{3}  = 'WBF_Z03_ZoneTempHSIMV';

inNms{1} = 'WBF_Z01_OilControl_FICHSIMV';
inNms{2} = 'WBF_Z01_CombAir_FICHSIMV';
inNms{3} = 'WBF_Z01_ZoneTempHSIWSP';
inNms{4} = 'WBF_Z02_OilControl_FICHSIMV';
inNms{5} = 'WBF_Z02_CombAir_FICHSIMV';
inNms{6} = 'WBF_Z02_ZoneTempHSIWSP';
inNms{7} = 'WBF_Z03_OilControl_FICHSIMV';
inNms{8} = 'WBF_Z03_CombAir_FICHSIMV';
inNms{9} = 'WBF_Z03_ZoneTempHSIWSP';
inNms{10} = 'SU_IML_GB6_SGNHSIValue';
inNms{11} = 'SU_UML_GB30_SGNHSIValue';
inNms{12} = 'WBF_MainExhaust_ExhaustFlow_FICHSIMV';
inNms{13} = 'WBF__PC027HSIMV';
inNms{14} = 'WBF_Z01_ColdCAirOutputHSIMV';
inNms{15} = 'WBF_Z02_ColdCAirOutputHSIMV';
inNms{16} = 'WBF_Z03_ColdCAirOutputHSIMV';
inNms{17} = 'WBF_Z01_O2_QICHSIMV';
inNms{18} = 'WBF_Z02_O2_QICHSIMV';
inNms{19} = 'WBF_Z03_O2_QICHSIMV';

% specify which inputs/outputs affect prediction of other inputs/outputs
% 0 means that there is no effect
% positive number n at (i,j) means that n past values of input/output 
% property j affect the prediction of i-th output
% dim( outCross ) = Ny x Ny
% dim( inCross ) = Ny x Nu

outCross = [4 2 0;...
            2 4 2;...
            0 2 4];
        
inCross = [ 3 3 0 0 0 0 0 0 0 3 0 0 0 3 0 0 3 0 0 ; ...
            0 0 0 3 3 0 0 0 0 0 0 0 0 0 3 0 0 3 0; ...
            0 0 0 0 0 0 3 3 0 0 3 0 0 0 0 3 0 0 3; 
          ];

        
%set the partition of the data used for training/testing
dataPartition.Training.start = 5000;
dataPartition.Training.end   = 30000;


dataPartition.Test.start = 75000;
dataPartition.Test.end = 85000;

%determine number of inputs and outputs
Ny = length(targetNms);
Nu = length(targetNms);

%various constants               
alpha = 0.1; %inptu to AVAR algorithm
K = 10; %how many steps ahead

%switch binary values for door openings (just for plots)
MEDFOS.SU_IML_GB6_SGNHSIValue = abs(MEDFOS.SU_IML_GB6_SGNHSIValue - 1);
MEDFOS.SU_UML_GB30_SGNHSIValue = abs(MEDFOS.SU_UML_GB30_SGNHSIValue - 1);

scaleHist=2; %scale the number of past values used in regression
polyStruct.inCross = scaleHist*inCross;
polyStruct.outCross = scaleHist*outCross;

data = medfosData( MEDFOS, targetNms, inNms, dataPartition, polyStruct);
data.scaleData();
data.createTrainingMat();


methodSelector = 2; %ordinary least squares

lassoMap = [];
if methodSelector==1
    [predictM, res] = data.createPredictModel('L1', 0);
elseif methodSelector==2
    [predictM, res] = data.createPredictModel('L2', 0);
elseif methodSelector==3
    [predictM, res] = data.createPredictModel('huber', c);
elseif methodSelector==4
    [predictM, res] = data.createPredictModel('L1S',  c );
elseif methodSelector==5
    [predictM, res] = data.createPredictModel('L2S', c);
else
    [predictM, res, lassoMap, sparsityInfo] = data.createPredictModelLasso( lambda );
end

[ Analytics, predictedTrajectory, measuredTrajectory]  = data.predictAhead( predictM, K, lassoMap );
perfIndice = data.calcPerformanceIndices( Analytics, alpha );



