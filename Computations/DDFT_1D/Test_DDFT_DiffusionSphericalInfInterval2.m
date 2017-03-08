function [output,optsNum,optsPhys,optsPlot] = Test_DDFT_DiffusionSphericalInfInterval2(doHI)

    if(nargin==0)
        doHI = true;
    end

    Phys_Area = struct('shape','InfSpectralLineSpherical','N',200,'L',4);
    Plot_Area = struct('N',200,'yMin',0,'yMax',12);
    Fex_Num   = struct('Fex','FMT','N',100);
    HI_Num    = struct('N',100,'L',2, ...
                       'HI11','noHISpherical',...
                       'HI12','RotnePrager12Spherical', ...
                       'HIPreprocess', 'RotnePragerPreprocessSpherical');
    
    tMax = 2.5;
    
    DDFTCode = 'DDFT_Diffusion_1D_Spherical';
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode',DDFTCode,...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    Vm = 0.1;

    ZS = [0.5;0.25];
    Z = repmat(ZS',Phys_Area.N,1);
    
    sigmaS  = [1 1;1 1];
    sigmaHS = [0.5 0.5; 0.5 0.5];
    
    V1       = struct('V1DV1','V1_Test_Spherical2','Vm',Vm,'Z',Z);
    V2       = struct('sigmaS',sigmaS);
        
    HI       = struct('sigmaHS',sigmaHS,'sigmaS',sigmaS);
    
    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',[25;25], ...
                           'mS',[1;1],'gammaS',[2;2]);
    
    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end                       
                       
    lineColourDDFT={{'m','b','r','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);  
end