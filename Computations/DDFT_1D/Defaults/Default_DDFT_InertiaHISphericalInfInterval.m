function [optsNum,optsPhys,optsPlot] = Default_DDFT_InertiaHISphericalInfInterval

    Phys_Area = struct('N',200,'L',4);
    Plot_Area = struct('N',200,'yMin',0,'yMax',8);
    Fex_Num   = struct('Fex','FMT','N',100);
    
    HI_Num    = struct('N',100,'L',4);
    
    tMax = 2.5;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'HINum',HI_Num,...
                         'DDFTCode','DDFT_InertiaSphericalInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);
    
    Vm = 0.1;

    ZS = [0.5;0.25];
    Z = repmat(ZS',Phys_Area.N,1);
    
    V1       = struct('V1DV1','APSHS','Vm',Vm,'Z',Z);

    sigmaS  = [1;1];
    sigmaHS = [0.5;0.5];
    
    HI       = struct('HI11','JeffreyOnishi11Spherical',...
              'HI12','JeffreyOnishi12Spherical', ...
              'HIPreprocess', 'JeffreyOnishiPreprocessSpherical', ...
              'sigma',[1 1;1 1],'sigmaH',[0.5 0.5; 0.5 0.5]); 
    
    optsPhys = struct('V1',V1,'V2',[],'HI',HI,'kBT',1,'nParticlesS',[25;25], ...
                           'mS',[1;1],'gammaS',[2;2],'sigmaS',sigmaS, ...
                           'sigmaHS',sigmaHS);
    
                      
                       
    lineColourDDFT={{'m','g','r','b'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
end