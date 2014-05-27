function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionHalfSpace_MF


    Phys_Area = struct('N',[30;30],'L1',2,'L2',2,'y2wall',0,'N2bound',10,'h',1,'L2_AD',1);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);
                   

    Fex_Num = struct('Fex','Meanfield','N',[20,20],'L',1);                   
                   
                   
    
    tMax = 20;

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','DDFT_DiffusionHalfSpace',...
                     'plotTimes',0:tMax/100:tMax);


    V1       = struct('V1DV1','V1_HalfSpace','V0',0.01,'V0Add',0.1,'sigma1',2,'sigma2',2,...
                        'tau',1,'y10',2,'y20',1);

    V2       = struct('V2DV2','Gaussian','alpha',1,'epsilon',0);
    
    optsPhys = struct('V1',V1,'V2',V2, ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',10); 
                  
    optsPlot.doDDFTPlots=true;
                  
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);                 

end                 

