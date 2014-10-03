function [optsNum,optsPhys,optsPlot,name] = Test_DDFT_FMT_InfCapillary()

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('N',[3;40],'L1',2,'y2Min',0.,'y2Max',3);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',4.5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','DDFT_FMT_InfCapillary',...
                     'plotTimes',0:7/100:7,...
                     'name','default');                                                  

    sigmaS = [1];%   1.1 ;
              %1.1 1.2 ];        
        
    V1       = struct('V1DV1','zeroPotential',...  %Vext_Cart_6
                      'V0',0.1,'V0r',1,'alphar',1,'tau',1,'rV',1);                                         
                  
    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 
    
    optsPhys = struct('V1',V1,'V2',V2,...                                            
                      'kBT',1,'eta',0.4257,...
                      'nParticlesS',10); %[10;10]
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    
    if(nargout == 0)
        f(optsPhys,optsNum,optsPlot);                 
    end

end                 
