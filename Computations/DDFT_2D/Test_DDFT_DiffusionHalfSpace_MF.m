function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionHalfSpace_MF


    Phys_Area = struct('shape','HalfSpace','N',[40;30],'L1',8,'L2',8,'y2wall',0);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0,'y2Max',5,'N2',100);
                   

    Fex_Num = struct('Fex','Meanfield','N',[20,20],'L',2,'L2',2);                   
                   
    tMax = 5;

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',Fex_Num,...                    
                     'plotTimes',0:tMax/100:tMax);


    V1       = struct('V1DV1','V1_HalfSpace','V0',0.2,'V0Add',1,'tau',1, ...
                      'sigma1',5,'sigma2',2);

    V2       = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
    
    optsPhys = struct('V1',V1,'V2',V2, ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',50); 
                  
    optsPlot.doDDFTPlots=true;
    
    config = v2struct(optsPhys,optsNum);
    
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium();
    EX.ComputeDynamics();    
    EX.PlotDynamics();
end                 

