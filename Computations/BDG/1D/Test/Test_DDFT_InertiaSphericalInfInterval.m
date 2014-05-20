function [output,optsPhys,optsNum,optsPlot] = Test_DDFT_InertiaSphericalInfInterval(doHI)

    if(nargin==0)
        doHI = true;
    end

    Phys_Area = struct('N',200,'L',4);
    Plot_Area = struct('N',200,'yMin',0,'yMax',8);
    Fex_Num   = struct('Fex','FMT','N',100);  
    HI_Num    = struct('N',100,'L',4, ...
                       'HI11','JeffreyOnishi11Spherical',...
                       'HI12','JeffreyOnishi12Spherical', ...
                       'HIPreprocess', 'JeffreyOnishiPreprocessSpherical');
    
    tMax = 2.5;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_InertiaSphericalInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);
    
    Vm = 0.1;

    ZS = [0.5;0.25];
    Z = repmat(ZS',Phys_Area.N,1);
    
    sigmaS  = [1 1;1 1];
    sigmaHS = [0.5 0.5; 0.5 0.5];

    V1       = struct('V1DV1','APSHS','Vm',Vm,'Z',Z);
    V2       = struct('sigmaS',sigmaS);
        
    HI       = struct('sigmaHS',sigmaHS,'sigmaS',sigmaS); 
    
    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',[25;25], ...
                           'mS',[1;1],'gammaS',[2;2]);
    
    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end                       
                       
    lineColourDDFT={{'m','g','r','b'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum,optsPlot);  
end