function output = TestFMT_DDFT_DiffusionInfSpace_HI(nexec)

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',... %'Rosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
    
    HI_Num    = struct('N',[20;20],'L',2);                   
                   
    %Sub_Area  = Phys_Area;
    
    tMax = 2;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'HINum',HI_Num,...
                     'DDFTCode','DDFT_DiffusionInfSpace_HI_NSpecies',...
                     'plotTimes',0:tMax/100:tMax);
    
    sigmaS = 1;
        
    V1       = struct('V1DV1','infHIDiffusion2',...
                      'V0',0.01,'V0add',2,'tau',0.1,'sigma1Add',5,'sigma2Add',10, ...
                      'y10',-3,'y20',0);                                         
                  
    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 
    
    HI       = struct('HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D', ...
                      'sigma',sigmaS,'sigmaH',sigmaS/2);
    
    optsPhys = struct('V1',V1,'V2',V2,'HI',HI,...                                            
                      'kBT',1,...
                      'nParticlesS',50); 
                 
    lineColourDDFT={{'b','r','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    if((nargin == 0) || ~nexec)
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        output = f(optsPhys,optsNum,optsPlot);                 
    end

end                 

