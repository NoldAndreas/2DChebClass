function [optsNum,optsPhys,optsPlot] = TestFMT_DDFT_DiffusionInfSpace_HI_NSpecies(nexec)

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N',[20,20],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',... %'Rosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
    
    HI_Num    = struct('N',[20;20],'L',2);                   
                   
    %Sub_Area  = Phys_Area;
    
    tMax = 10;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'HINum',HI_Num,...
                     'DDFTCode','DDFT_DiffusionInfSpace_HI_NSpecies',...
                     'plotTimes',0:tMax/100:tMax);

	sigma1 = 1;
    sigma2 = 2;
    sigma12 = (sigma1+sigma2)/2;
                 
    sigmaS = [ sigma1   sigma12 ;
               sigma12  sigma2 ];        
        
    V1       = struct('V1DV1','rotating2cart',...
                      'V0',0.01,'V0r',1,'alphar',3,'tau',0.5,'rV',1);                                         
                  
    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 
    
    HI       = struct('HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D', ...
                      'sigma',sigmaS,'sigmaH',sigmaS/2);
    
    optsPhys = struct('V1',V1,'V2',V2,'HI',HI,...                                            
                      'kBT',1,...
                      'nParticlesS',[10;10]); 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    if((nargin == 0) || ~nexec)
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        f(optsPhys,optsNum,optsPlot);                 
    end

end                 

