function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionHalfSpace_FMT_Triangle(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','HalfSpace_FMT','N',[20;20],'L1',2,'L2',2, ...
                       'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90);                    
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
                   
    HI_Num    = struct('N',[10;10],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');  
    
    tMax = 0.25;

    DDFTCode = 'DDFT_DiffusionHalfSpace_FMT_skewed';
%     DDFTCode = 'DDFT_Diffusion_2D';    
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode',DDFTCode,...
                     'plotTimes',0:tMax/100:tMax, ...
                     'Accuracy_Averaging',1e-6);

    sigmaS  = 1;
    sigmaHS = 0.5;
    
    V1       = struct('V1DV1','V1_Triangle',...
                      'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
                      'y10',-1,'y20',1.5,'y11',1,'y21',2,'y12',0,'y22',2.5); 

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

