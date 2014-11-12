function [optsNum,optsPhys] = DDFT_DiffusionBox_2Species_1()

 
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
                     'plotTimes',0:1:75);
                 
	epsilonS=  [ 1   -0.1 ;
                 -0.1  1 ] ;
             
    M = 900;

    epsilon_w1     = repmat([1,1],M,1);
    epsilon_w1_end = repmat([1,0],M,1);
    epsilon_w2     = repmat([0,0],M,1);
    epsilon_w2_end = repmat([0,1],M,1);                
    
    V1        = struct('V1DV1','Vext_Cart_2Species_1',...
                      'V0',0.0,'y10',4,'y20',4,'tau',5,...
                      'epsilon_w1',epsilon_w1,'epsilon_w1_end',epsilon_w1_end,...
                      'epsilon_w2',epsilon_w2,'epsilon_w2_end',epsilon_w2_end);    
    V2       = struct('V2DV2','Phi2DLongRange','epsilon',epsilonS);                                     
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,...
                     'HSBulk','CarnahanStarling',...                     
                     'nParticlesS',[30;30]);


    config = v2struct(optsPhys,optsNum);
    %****************************
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium(0.5);%EX.optsPhys.rhoLiq_sat);
    EX.ComputeDynamics();           
    EX.PlotDynamics();    
    
    

end                 

