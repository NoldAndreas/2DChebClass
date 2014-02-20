                 
function [optsNum,optsPhys] = Default1_DDFT_DiffusionBox()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',10,'N',[25;25],...
                       'y2Min',0,'y2Max',10);
    Plot_Area = Phys_Area; 
    Plot_Area.N1 = 100; Plot_Area.N2 = 100;
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',1,'y1Max',3,'N1',20,...
                      'y2Min',1,'y2Max',2,'N2',20);
    
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.2:3,...
                     'DDFTCode','DDFT_DiffusionBox');    
                      
    V1 = struct('V1DV1','Vext_Cart_1','V0',0.1,'grav',2,'y10',2,'y20',2);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
    
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'nParticlesS',50);

    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);    

end                 

