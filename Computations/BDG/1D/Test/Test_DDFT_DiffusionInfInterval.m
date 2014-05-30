function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionInfInterval(doHI)

    if(nargin==0)
        doHI = true;
    end

    Phys_Area = struct('shape','InfSpectralLine','N',200,'L',30);
    Plot_Area = struct('N',200,'yMin',-20,'yMax',20);
    Fex_Num   = struct('Fex','Percus','N',100);
    
    HI_Num    = struct('N',100,'L',2, ...
                      'HI11','noHI','HI12','RotnePrager12', ...
                      'HIPreprocess', 'RotnePragerPreprocess');
    
    tMax = 25;
    
    %DDFTCode = 'DDFT_DiffusionInfInterval';
    DDFTCode = 'DDFT_Diffusion_1D_Planar';
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode',DDFTCode,...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    alpha0 = 0.01;
    alphaT = 15;
    beta   = 40;
     
    V1       = struct('V1DV1','V1_Test_Infinite','alpha0',alpha0,'alphaT',alphaT,'beta',beta);
    
    V2       = struct('sigmaS',1);
    
    HI       = struct('sigmaHS',0.5);

    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',8, ...
                           'mS',1,'gammaS',2);
    
    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
        
    lineColourDDFT={{'g','r','b'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);  
    
end