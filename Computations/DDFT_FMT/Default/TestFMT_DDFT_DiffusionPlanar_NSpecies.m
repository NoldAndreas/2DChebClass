function [optsNum,optsPhys,optsPlot] = TestFMT_DDFT_DiffusionPlanar_NSpecies(nexec)

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N',[20,20],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',... %'Rosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','DDFT_DiffusionPlanar_NSpecies',...
                     'plotTimes',0:4/100:4);
                     
    sigmaS = [ 1   1.1 ;
              1.1 1.2 ];        
        
    V1       = struct('V1DV1','rotating2cart',...
                      'V0',0.01,'V0r',1,'alphar',2,'tau',1,'rV',1);                                         
                  
    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 
    
    optsPhys = struct('V1',V1,'V2',V2,...                                            
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

