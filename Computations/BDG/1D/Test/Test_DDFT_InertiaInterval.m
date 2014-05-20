function [output,optsNum,optsPhys,optsPlot] = Test_DDFT_InertiaInterval()

    yMax = 13.29;

    Phys_Area = struct('N',100,'yMin',0,'yMax',yMax);
    Plot_Area = struct('N',200,'yMin',0,'yMax',yMax);
    Fex_Num   = struct('Fex','Percus','N',100);
    
    tMax = 10;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_InertiaInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax,...
                         'doPlots',true);

    alpha = repmat([5,5],Phys_Area.N,1);
    beta = repmat([0.5,0.5],Phys_Area.N,1);
    y0 = repmat([4,8],Phys_Area.N,1);
    
    sigmaS  = [1 1;1 1];
    
    V1       = struct('V1DV1','bump1D','alpha',alpha,'beta',beta,'y0',y0);
    V2       = struct('sigmaS',sigmaS);      
     
    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',[5;5], ...
        'mS',[1;1],'gammaS',[5;5]);
    
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);  
    
end