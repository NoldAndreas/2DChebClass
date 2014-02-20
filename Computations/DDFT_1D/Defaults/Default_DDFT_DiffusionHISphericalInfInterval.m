function [optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionHISphericalInfInterval

    Phys_Area = struct('N',200,'L',4);
    Plot_Area = struct('N',200,'yMin',0,'yMax',12);
    Fex_Num   = struct('Fex','FMT','N',100);
    HI_Num    = struct('N',100,'L',2);
    
    tMax = 2.5;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'HINum',HI_Num,...
                         'DDFTCode','DDFT_DiffusionSphericalInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    Vm = 0.1;

    ZS = [0.5;0.25];
    Z = repmat(ZS',Phys_Area.N,1);
    
    V1       = struct('V1DV1','APSHS','Vm',Vm,'Z',Z);

    HI       = struct('HI11','noHISpherical',...
                      'HI12','RotnePrager12Spherical', ...
                      'HIPreprocess', 'RotnePragerPreprocessSpherical', ...
                      'sigma',[1 1;1 1],'sigmaH',[0.5 0.5;0.5 0.5]);
    
    optsPhys = struct('V1',V1,'V2',[],'HI',HI,'kBT',1,'nParticlesS',[25;25], ...
                           'mS',[1;1],'gammaS',[2;2],'sigmaS',[1;1],...
                           'sigmaHS',[0.5;0.5]);
    
    lineColourDDFT={{'m','b','r','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
end