function [optsNum,optsPhys,optsPlot] = TestFMT_DDFT_DiffusionPlanar_1Species(doHI)

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N',[20;20],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    xPlotLim = 10;
    yPlotLim = 10;
    
    Plot_Area = struct('y1Min',-xPlotLim,'y1Max',xPlotLim,'N1',100,...
                       'y2Min',-yPlotLim,'y2Max',yPlotLim,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',... %'Rosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);
    
    tMax = 0.7;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','DDFT_DiffusionPlanar_NSpecies',...
                     'plotTimes',0:tMax/100:tMax);
                        
    sigmaS = 1;        
        
%     V1       = struct('V1DV1','quadQuad',...
%                       'V0x',0.01,'V0y',0.02,'tau',0.01,'V0Add',3,'A0',1,'B0x',3,'B0y',NaN); 

    V1       = struct('V1DV1','quadQuad',...
                      'V0x',0.01,'V0y',0.02,'tau',0.01,'V0Add',3,'A0',1,'B0x',4,'B0y',NaN); 


    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 
    
    
    
    optsPhys = struct('V1',V1,'V2',V2,...                                            
                      'kBT',1,...
                      'nParticlesS',10); 
   
    optsPhys.doHI = doHI;
                  
    if(doHI)
        HIShapeParams = struct('RMinS',1,'N',[20;20],'L',2);
        optsNum.HIShapeParams = HIShapeParams;
        HI       = struct('sigmaHS',1);
        optsPhys.HI=HI;
    end
                  
                  
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
%    if((nargin == 0) || ~nexec)
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        f(optsPhys,optsNum,optsPlot);                 
%    end

end                 

