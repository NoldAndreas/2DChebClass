function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionHalfSpace_FMT(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('N',[30;30],'L1',2,'L2',2,'y2wall',0,'N2bound',10,'h',1,'L2_AD',1);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
                   
    HI_Num    = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');  
    
    tMax = 10;

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','DDFT_DiffusionHalfSpace_FMT',...
                     'plotTimes',0:tMax/100:tMax, ...
                     'Accuracy_Averaging',1e-6);

    sigmaS  = 1;
    sigmaHS = 0.5;

    V1       = struct('V1DV1','V1_HalfSpace_FMT','V0',0.1,'V0Add',1,'sigma1',1,'sigma2',1,...
                        'tau',1,'y10',2,'y20',2);

    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);   

    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1,'V2',V2, ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',10); 

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;
                  
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);                 

end                 

