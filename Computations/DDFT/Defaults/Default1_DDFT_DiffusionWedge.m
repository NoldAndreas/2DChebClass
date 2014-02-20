function [optsNum,optsPhys] = Default1_DDFT_DiffusionWedge()
    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',3,'L1',2,'N',[30,20],...
                       'y2Min',-pi*2/3,'y2Max',pi*2/3);
    Plot_Area       = Phys_Area; 
    Plot_Area.N1    = 100; 
    Plot_Area.N2    = 100; 
    Plot_Area.y1Max = 3;
    
%    Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',1,'y1Max',3,'N1',20,...
                      'y2Min',-pi*1/3,'y2Max',pi*1/3,'N2',20);
    
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:2,...
                     'DDFTCode','DDFT_DiffusionWedge');
                 
    V1 = struct('V1DV1','Vext_Pol_2','V0',0.1,'grav',1,'r0',1);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'nParticlesS',50);

    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);   
end