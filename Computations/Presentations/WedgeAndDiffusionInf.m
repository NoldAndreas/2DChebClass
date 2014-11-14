%To Show

y1Min = 3;
y1Max = 6;
y2Min = -1;
y2Max = 5;

N  = [30;30];

st = v2struct(y1Min,y1Max,y2Min,y2Max,N);

BX = Box(st);

stPlot = st;
stPlot.N = [100;100];

BX.ComputeAll(stPlot);

y = cos(BX.Pts.y1_kv).*sin(BX.Pts.y2_kv);

subplot(2,2,1);  BX.plot(y);
subplot(2,2,2);  BX.plot(y,'contour');
subplot(2,2,3);  BX.plot(BX.Diff.Dy2*y);
subplot(2,2,4);  BX.plot(BX.Diff.Lap*y);

SimulationWedge();
DiffusionAdvectionPolarInfinity();

%Test0_DDFT_DiffusionBox_2Species();