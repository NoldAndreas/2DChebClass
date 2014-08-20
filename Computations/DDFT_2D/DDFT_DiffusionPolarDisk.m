function [optsNum,optsPhys] = DDFT_DiffusionPolarDisk()        

    Phys_Area = struct('shape','Disc','N',[20,20],'R',2);

    Plot_Area = struct('y1Min',0,'y1Max',2,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
                   
	V2Num = struct('Fex','Meanfield','N',[20,20]);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',V2Num,...
                     'plotTimes',0:0.2:4);                            

    V1 = struct('V1DV1','Vext_Cart_3','V0',0.5,'grav',1);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'nParticlesS',50);  
    
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                 
    AddPaths();
    DDFTDynamics(optsPhys,optsNum,optsPlot);
   
end
    