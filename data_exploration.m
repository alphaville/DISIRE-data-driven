%% Intro

load data/medfos.mat
set(0,'DefaultAxesFontSize',16)
%%
%stem(find(x~=0), x(x~=0).*(1-0.15*randn(5,1)),'m','LineWidth',2);
%hold on
%stem(find(x~=0), x(x~=0),'s','LineWidth',2);
%hold on



%% Temperatures: T1-T2
dat=[MEDFOS.WBF_Z01_ZoneTempHSIMV,MEDFOS.WBF_Z02_ZoneTempHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(T_1, T_2)');
grid on
view(3);

savefig('plots/data_expl_T12.fig');
print('-painters','-deps','plots/data_expl_T12.epsc')

%% Temperatures: T1-T3
dat=[MEDFOS.WBF_Z01_ZoneTempHSIMV,MEDFOS.WBF_Z03_ZoneTempHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(T_1, T_3)');
grid on
view(3);

savefig('plots/data_expl_T13.fig');
print('-painters','-deps','plots/data_expl_T13.epsc')


%% Temperatures: T2-T3
dat=[MEDFOS.WBF_Z02_ZoneTempHSIMV,MEDFOS.WBF_Z03_ZoneTempHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(T_2, T_3)');
grid on
view(3);

savefig('plots/data_expl_T23.fig');
print('-painters','-deps','plots/data_expl_T23.epsc')

%% AIR: AIR1-AIR2
dat=[MEDFOS.WBF_Z01_CombAir_FICHSIMV,MEDFOS.WBF_Z02_CombAir_FICHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(AIR_1, AIR_2)');
grid on
view(3);

savefig('plots/data_expl_air12.fig');
print('-painters','-deps','plots/data_expl_air12.epsc')

%% AIR: AIR1-AIR3
dat=[MEDFOS.WBF_Z01_CombAir_FICHSIMV,MEDFOS.WBF_Z03_CombAir_FICHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(AIR_1, AIR_3)');
grid on
view(3);
savefig('plots/data_expl_air13.fig');
print('-painters','-deps','plots/data_expl_air13.epsc')

%% AIR: AIR2-AIR3
dat=[MEDFOS.WBF_Z01_CombAir_FICHSIMV,MEDFOS.WBF_Z03_CombAir_FICHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(AIR_2, AIR_3)');
grid on
view(3);

savefig('plots/data_expl_air23.fig');
print('-painters','-deps','plots/data_expl_air23.epsc')

%% Pressure - Temperature 1

dat=[MEDFOS.WBF__PC027HSIMV,MEDFOS.WBF_Z03_CombAir_FICHSIMV];
hist3(dat);
hold on
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);

h = pcolor(xb,yb,n1/100);

h.ZData = ones(size(n1)) * -max(max(n));


%colormap(hot) % heat map
title('(P, AIR_3)');
grid on
view(3);


savefig('plots/data_expl_PA3.fig');
print('-painters','-deps','plots/data_expl_PA3.epsc');


%% Tracking error
varx = MEDFOS.WBF_Z03_CombAir_FICHSIMV;
varxsp = MEDFOS.WBF_Z03_CombAir_FICHSIWSP;
time=(1:length(varx))./(360);
plot(time,varx,'LineWidth',2)
hold on
plot(time,varxsp)
grid on
axis([0, max(time), min(varx), 1.1*max(varx)])
xlabel('time [hr]');
ylabel('Air_3');
savefig('plots/Air3.fig');
print('-painters','-deps','plots/Air3.epsc');


%% Exhaust flow

plot(time, MEDFOS.WBF__PC027HSIMV)
axis tight
grid on
xlabel('time [hr]');
ylabel('Pressure')
savefig('plots/pressure.fig');
print('-painters','-deps','plots/pressure.epsc')

%% Doors 

door1 = MEDFOS.SU_IML_GB6_SGNHSIValue;
door2 = MEDFOS.SU_UML_GB30_SGNHSIValue;
time=(1:length(door1))./(360);
subplot(511)
plot(time,~door1,'LineWidth',1)
axis tight
ylabel('Door A');
subplot(512)
plot(time,~door2,'LineWidth',1)
ylabel('Door B');
axis tight
subplot(513)
plot(time, 100*(MEDFOS.WBF_Z01_CombAir_FICHSIMV-MEDFOS.WBF_Z01_CombAir_FICHSIWSP)./ MEDFOS.WBF_Z01_CombAir_FICHSIWSP);
ylabel('Air_1');
axis tight
subplot(514)
plot(time, 100*(MEDFOS.WBF_Z02_CombAir_FICHSIMV-MEDFOS.WBF_Z02_CombAir_FICHSIWSP)./MEDFOS.WBF_Z02_CombAir_FICHSIWSP);
ylabel('Air_2');
axis tight
subplot(515)
plot(time, 100*(MEDFOS.WBF_Z03_CombAir_FICHSIMV-MEDFOS.WBF_Z03_CombAir_FICHSIWSP)./MEDFOS.WBF_Z03_CombAir_FICHSIWSP);
ylabel('Air_3');
axis tight
reply = input('strike any key to proceed:','s');
savefig('plots/doors_effect_air.fig');
print('-painters','-deps','plots/doors_effect_air.epsc');
%%

hold off
subplot(511)
plot(time,~door1,'LineWidth',1)
axis tight
ylabel('Door A');
subplot(512)
plot(time,~door2,'LineWidth',1)
ylabel('Door B');
axis tight
subplot(513)
plot(time, 100*(MEDFOS.WBF_Z01_ZoneTempHSIMV-MEDFOS.WBF_Z01_ZoneTempHSIWSP)./ MEDFOS.WBF_Z01_ZoneTempHSIWSP);
ylabel('T_1');
axis tight
subplot(514)
plot(time, 100*(MEDFOS.WBF_Z02_ZoneTempHSIMV-MEDFOS.WBF_Z02_ZoneTempHSIWSP)./MEDFOS.WBF_Z02_ZoneTempHSIWSP);
axis tight
ylabel('T_2');
subplot(515)
plot(time, 100*(MEDFOS.WBF_Z03_ZoneTempHSIMV-MEDFOS.WBF_Z03_ZoneTempHSIWSP)./MEDFOS.WBF_Z03_ZoneTempHSIWSP);
axis tight
ylabel('T_3');
savefig('plots/doors_effect_temp.fig');
print('-painters','-deps','plots/doors_effect_temp.epsc');


%% Fuel
time=(1:length(varx))./(360);
plot(time, MEDFOS.WBF_Z03_OilControl_FICHSIMV)
xlabel('time [hr]');
ylabel('Fuel #3')
axis tight
grid on
savefig('plots/fuel3.fig');
print('-painters','-deps','plots/fuel3.epsc');


