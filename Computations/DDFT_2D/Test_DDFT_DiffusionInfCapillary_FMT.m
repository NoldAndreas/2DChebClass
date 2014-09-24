function EX = Test_DDFT_DiffusionInfCapillary_FMT(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','InfCapillary_FMT','N',[30;30],'L1',2,'L2',2, ...
                       'y2Min',0,'y2Max',10,'N2bound',10,'L2_AD',1);
    
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
                     'plotTimes',0:tMax/100:tMax);

    sigmaS  = 1;
    sigmaHS = 0.5;

    V1       = struct('V1DV1','V1_HalfSpace_FMT','V0',0.1,'V0Add',1,'sigma1',1,'sigma2',1,...
                        'tau',1,'y10',2,'y20',2);    

    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1, ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',10,'sigmaS',sigmaS); 

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;
    
    EX = DDFTDynamics(optsPhys,optsNum,optsPlot);
    
end                 
