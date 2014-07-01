                 
function [optsNum,optsPhys] = DDFT_DiffusionBox_2Phase_1()

    %Numerical Parameters    
    Phys_Area = struct('shape','Box','N',[25,25],...
                       'y1Min',0,'y1Max',10,...
                       'y2Min',-6,'y2Max',6);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',-6,'y2Max',6,'L2',10,'N2',100);
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('shape','Box','y1Min',5,...
                      'y1Max',10,'N',[20,20],...
                      'y2Min',0,'y2Max',5);
                  
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'V2Num',FexNum,...
                     'plotTimes',0:1:20);                    
                     
    V1 = struct('V1DV1','Vext_Cart_5',...
                      'V0',0.02,'epsilon_w',1,'epsilon_w_end',0,...
                      'y10',2,'y20',0,'tau',1);
                  
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);

    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,'nParticlesS',40);                                  
               
                  
	config = v2struct(optsPhys,optsNum);
        
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium( EX.optsPhys.rhoGas_sat);
    EX.ComputeDynamics();    

end                 
