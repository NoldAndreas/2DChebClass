clear all
close all

Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);

Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                    'y2Min',-5,'y2Max',5,'N2',100);
                   
IDC = InfSpace_FMT(Phys_Area);

[Pts,Diff,Int,Ind,~]  = IDC.ComputeAll(Plot_Area);

y1_0 = 0;
y2_0 = 0;

optsPhys.sigmaHS = 1;

HI0 = RP12_2D(Pts.y1_kv-y1_0,Pts.y2_kv-y2_0,optsPhys);

HI0 = reshape(HI0,[size(HI0,1),4]);

IDC.doPlots(HI0);

optsPhys.HIfn = str2func('RP12_2D');
optsPhys.R    = optsPhys.sigmaHS/2;

HIm = MollifyHIFunction(Pts.y1_kv-y1_0,Pts.y2_kv-y2_0,optsPhys);

HIm = reshape(HIm,[size(HIm,1),4]);

figure
IDC.doPlots(HIm);

HIdiff = HI0 - HIm;

figure
IDC.doPlots(HIdiff);