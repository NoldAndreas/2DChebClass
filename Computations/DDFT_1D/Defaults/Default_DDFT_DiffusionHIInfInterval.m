function [optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionHIInfInterval

    Phys_Area = struct('N',200,'L',30);
    Plot_Area = struct('N',200,'yMin',-20,'yMax',20);
    Fex_Num   = struct('Fex','Percus','N',100);
    
    HI_Num    = struct('N',100,'L',2);
    
    tMax = 25;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'HINum',HI_Num,...
                         'DDFTCode','DDFT_DiffusionHIInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    alpha0 = 0.01;
    alphaT = 15;
    beta   = 40;
     
    V1       = struct('V1DV1','oscBump1D','alpha0',alpha0,'alphaT',alphaT,'beta',beta);
    
    HI       = struct('HI11','noHI','HI12','RotnePrager12', ...
                      'HIPreprocess', 'RotnePragerPreprocess', ...
                      'sigma',1,'sigmaH',0.5);

    
    optsPhys = struct('V1',V1,'V2',[],'HI',HI,'kBT',1,'nParticlesS',8,'sigmaS',1, ...
                           'mS',1,'gammaS',2);
    
    lineColourDDFT={{'g','r','b'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
end