set(0,'DefaultAxesFontSize',16)

% Length to plot
t0 = 1;
N = 111773;

% Establish a connection to the server
db = DisireConnection.instance();

% The three combustion air flows

idx1=db.getFineData(t0,N,'WBF_Z01_CombAir_FIC:HSI.MV');
idx2=db.getFineData(t0,N,'WBF_Z02_CombAir_FIC:HSI.MV');
idx3=db.getFineData(t0,N,'WBF_Z03_CombAir_FIC:HSI.MV');

% Total flow (F1+F2+F3)
total_flow = idx1.value+idx2.value+idx3.value;

% exhaust
idxe=db.getFineData(t0,N,'WBF_MainExhaust_ExhaustFlow_FIC:HSI.MV');

% pressure
idxP=db.getFineData(t0,N,'WBF__PC027:HSI.MV');
%% PLOT
subplot(211);
plot(idx1.value); hold on;
plot(idx2.value);
plot(idx3.value);
plot(total_flow,'--','Linewidth',2);
plot(idxe.value,'Linewidth',3);
legend('F_1','F_2','F_3','F_{tot}','F_e');
ylabel('Air Flow');

subplot(212);
plot(idxP.value);
legend('P');
ylabel('Pressure');
xlabel('Time');