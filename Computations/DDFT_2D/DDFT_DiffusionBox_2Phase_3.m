function [config,res] = DDFT_DiffusionBox_2Phase_3()
    
    Phys_Area    = struct('shape','Box','N',[30,30],...
                       'y1Min',0,'y1Max',10,...
                       'y2Min',0,'y2Max',10);
    
    Plot_Area    = Phys_Area; 
    Plot_Area.N1 = 100; 
    Plot_Area.N2 = 100;
        
    Sub_Area     = struct('shape','Box','N',[35,35],...
                      'y1Min',5,'y1Max',10,...
                      'y2Min',0,'y2Max',5);    
                  
    FexNum       = struct('Fex','Meanfield','N',[20,20],'L',1);                  
        
    optsNum      = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'V2Num',FexNum,...                     
                     'plotTimes',0:0.2:5);                            
                 
    V1 = struct('V1DV1','Vext_Cart_2Species_1',...
                     'V0',0.02,'grav',2,'y10',5,'y20',5,'tau',1,...                      
                     'epsilon_w1',1,'epsilon_w1_end',0,...
                     'epsilon_w2',0,'epsilon_w2_end',1);
                 
    V2 = struct('V2DV2','Phi2DLongRange','alpha',2,'epsilon',1);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,...
                     'HSBulk','CarnahanStarling',...                     
                     'nParticlesS',40);


    config = v2struct(optsPhys,optsNum);
    %****************************
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium([],struct('Iterative',true));
    % EX.ComputeEquilibrium(EX.optsPhys.rhoLiq_sat);
    EX.ComputeDynamics();           
    res.fig_handles = EX.PlotDynamics();                 
end                 
