function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionInfCapillary()

    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'L1',3,'N',[20;20],...
                       'y2Min',-2,'y2Max',2);
                   
	Plot_Area = struct('y1Min',-5,'y1Max',5,'y2Min',-2,'y2Max',2, ...
                       'N1',100,'N2',100);
    
    Sub_Area        = struct('y1Min',-3,'y1Max',3,'N',[20;20],...
                             'y2Min',-2,'y2Max',2);
    
    FexNum = struct('Fex','Meanfield','N',[20,20],'L1',6,'L2',1);
    %Conv      = struct('L1',6,'L2',1,'N',[30,30]);

    tMax = 5;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'FexNum',FexNum,...
                     'DDFTCode','DDFT_DiffusionInfCapillary', ...
                     'plotTimes',0:tMax/100:tMax);
                 
    V1 = struct('V1DV1','V1_InfCapillary','V0',0.1,'grav',1.0);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',1,'nParticlesS',50,'mS',1,'gammaS',1);
                 
    optsPlot.doDDFTPlots=true;

    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum);    

end       