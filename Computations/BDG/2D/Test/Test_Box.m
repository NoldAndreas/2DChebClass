function [optsNum,optsPhys] = Test_Box()

    Phys_Area = struct('N',[20;20], ...
                       'y1Min',0,'y1Max',10,...
                       'y2Min',0,'y2Max',10);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',0,'y2Max',10,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    
    tMax = 10;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',FexNum,...
                     'DDFTCode','DDFT_DiffusionBox',...
                     'plotTimes',0:tMax/100:tMax);
                 
                      
    V1 = struct('V1DV1','V1_Droplet','V0',1,'sigma',3,'tau',1,'y10',3,'y20',0);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.05);
    %V2 = struct('V2DV2','Phi2DLongRange','epsilon',0.1);
    
    
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',1,'mS',1,'gammaS',1,'nParticlesS',50);

    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);    

end                 

