function [optsNum,optsPhys] = Default1_DDFT_DiffusionPolarDisk()        


    Phys_Area = struct('y1Min',0,'y1Max',2,'L1',1,'N1',20,...
                       'y2Min',0,'y2Max',2*pi,'N2',20);

    Plot_Area = struct('y1Min',0,'y1Max',2,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'plotTimes',0:0.2:4,...
                     'DDFTCode','DDFT_DiffusionPolar');                            

    V1 = struct('V1DV1','Vext_Pol_3','V0',0.5,'grav',1);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'nParticlesS',50);  

    if(nargout == 0)
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        f(optsPhys,optsNum);                                             
    end
end
    