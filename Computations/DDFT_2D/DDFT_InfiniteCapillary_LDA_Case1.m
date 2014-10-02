function [optsNum,optsPhys] = DDFT_InfiniteCapillary_LDA_Case1()
    %Numerical Parameters    
    Phys_Area = struct('shape','InfCapillary',...
                       'y1Min',-inf,'y1Max',inf,'L1',3,'N',[20;20],...
                       'y2Min',-2,'y2Max',2);

	Plot_Area       = Phys_Area;
    Plot_Area.N1    = 100; Plot_Area.N2 = 100;
    Plot_Area.y1Min = -5;
    Plot_Area.y1Max = 5;
    
    Sub_Area        = struct('shape','Box',...
                            'y1Min',-3,'y1Max',3,'N',[20;20],...
                             'y2Min',-2,'y2Max',2);
        
    V2Num   = struct('Fex','Meanfield','L',1,'L1',6,'L2',1,'N',[30,30]);            
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.2:5,...
                     'V2Num',V2Num); 
                 
    V1 = struct('V1DV1','Vext_Cart_2','V0',0.1,'grav',1.0);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',-0.06);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'nParticlesS',50);
    
	DDFTDynamics(optsPhys,optsNum);

end       