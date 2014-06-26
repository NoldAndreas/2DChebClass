function [optsNum,optsPhys] = DDFT_DiffusionBox_2Phase_2()

    %Numerical Parameters    
    Phys_Area = struct('shape','Box','y1Min',0,'y1Max',10,'N',[30,30],...
                       'y2Min',0,'y2Max',10);
    
    Plot_Area = Phys_Area; 
    Plot_Area.N1 = 100; Plot_Area.N2 = 100;
        
    Sub_Area = struct('shape','Box','y1Min',5,'y1Max',10,'N',[35,35],...
                      'y2Min',0,'y2Max',5);    
                  
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);                  
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'V2Num',FexNum,...                     
                     'plotTimes',0:0.2:3);                            

    V1 = struct('V1DV1','Vext_Cart_3','V0',0.02,'grav',2,'y10',5,'y20',5);
    V2 = struct('V2DV2','Phi2DLongRange','alpha',2,'epsilon',1.5);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,...
                     'HSBulk','CarnahanStarling',...                     
                     'nParticlesS',40,'gammaS',1);
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;                 

    config = v2struct(optsPhys,optsNum);
    %****************************
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium();
    EX.ComputeDynamics();               
                 
end                 

