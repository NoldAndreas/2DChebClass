function output = Default_DDFT_DiffusionInfSpace_HS_HI

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);

    HI_Num    = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');                   
                   
    %Sub_Area  = Phys_Area;
    
    tMax = 0.5;

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'HINum',HI_Num,...
                     'DDFTCode','DDFT_DiffusionInfSpace',...
                     'plotTimes',0:tMax/100:tMax);

    sigmaS  = 1;
    sigmaHS = 0.5;                   

    V1       = struct('V1DV1','triangleDiffusion',...
                      'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
                      'y10',-1,'y20',-1,'y11',1,'y21',-1,'y12',0,'y22',0.5); 
                  
    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);   
    
    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1,'V2',V2,'HI',HI, ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',40); 
                 
    optsPlot.doDDFTPlots=true;
                  
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);                 

end                 

