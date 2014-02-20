function [optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionSphericalInfInterval

    Phys_Area = struct('N',200,'L',4);
    Plot_Area = struct('N',200,'yMin',0,'yMax',10);
    Fex_Num   = struct('Fex','Meanfield','N',100,'L',2);
    
    tMax = 12;
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_DiffusionSphericalInfInterval',...
                         'Tmax',tMax,'plotTimes',0:tMax/100:tMax, ...
                         'doPlots',true);

    Vm = 0.01;
    tSwitch = tMax;

    V1       = struct('V1DV1','APSG','Vm',Vm,'tSwitch',tSwitch);

    epsilon = [0.5,0.5;0.5,0.5];
    
    alpha1  = 0.2;
    alpha2  = 1;
    alpha12 = (alpha1+alpha2)/2;
    
    alpha   = [alpha1, alpha12; alpha12, alpha2]; 
    
    V2       = struct('V2DV2','GaussianSpherical','alpha',alpha,'epsilon',epsilon);   
     
    optsPhys = struct('V1',V1,'V2',V2,'kBT',1,'nParticlesS',[25;25], ...
                           'mS',[1;1],'gammaS',[2;2]);
    
    lineColourDDFT={{'r','b','m','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
end