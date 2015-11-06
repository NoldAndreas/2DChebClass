function [output,optsNum,optsPhys,optsPlot] = Test_NS_InfInterval(tMax)
    
    Phys_Area = struct('N',200,'L',30);
    Plot_Area = struct('N',200,'yMin',-20,'yMax',20);
    Fex_Num   = struct('Fex','Zero','N',100);
    %Fex_Num   = struct('Fex','Percus','N',100);
    %Fex_Num   = struct('Fex','Meanfield','N',100,'L',2);
    
    if(nargin==0)
        tMax = 25;
    end
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','NS_InfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    %alpha = repmat([3,3],Phys_Area.N,1);
    
    alpha0 = 0.01;
    alphaT = 15;
    beta   = 40;
    
    V1       = struct('V1DV1','V1_Test_NS','alpha0',alpha0,'alphaT',alphaT,'beta',beta);
      
    V2       = struct('V2DV2','zeroInteraction1D');   
     
    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',8, ...
                           'mS',1,'gammaS',2,'eta',0.1);
    
    lineColourDDFT={{'k','r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);  
end