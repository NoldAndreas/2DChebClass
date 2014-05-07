function [optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionHSSphericalInfInterval2

    Phys_Area = struct('N',200,'L',15);
    Plot_Area = struct('N',200,'yMin',0,'yMax',15);
    Fex_Num   = struct('Fex','FMT','N',100);
    
    tMax = 2.5;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_DiffusionSphericalInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);
    Vm = 0.1;

    ZS = [1;0.5];
    Z = repmat(ZS',Phys_Area.N,1);
    
    sigmaS  = [1 1;1 1];
    
    V1       = struct('V1DV1','APSHS','Vm',Vm,'Z',Z,'yI',2,'yt',8,'yt2',4);
    V2       = struct('sigmaS',sigmaS);
    
    optsPhys = struct('V1',V1,'V2',[],'kBT',1,'nParticlesS',[200;200], ...
                           'mS',[1;1],'gammaS',[10;10]);
    
    lineColourDDFT={{'m','b','r','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
end