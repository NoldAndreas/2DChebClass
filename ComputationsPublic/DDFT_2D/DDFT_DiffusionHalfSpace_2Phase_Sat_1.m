function [config,res] = DDFT_DiffusionHalfSpace_2Phase_Sat_1()

    %Numerical Parameters    
    Phys_Area = struct('shape','HalfSpace','N',[10,40],...
                       'y2Min',0,'L1',4,'L2',5);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,'N2',100,...
                       'y2Min',0,'y2Max',10);
        
    Sub_Area = struct('shape','Box','y1Min',-1,'y1Max',1,'N',[20,20],...
                      'y2Min',0,'y2Max',2);
                  
	V2Num    = struct('Fex','Meanfield','L',2,'L2',1.,'N',[40,40]);
        
    optsNum  = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:6,...
                     'V2Num',V2Num);
                     %'DDFTCode','DDFT_DiffusionHalfSpace_2Phase_Sat');
                      
 
    V1 = struct('V1DV1','Vext_Cart_5',...
                      'V0',0.0,'y10',2,'y20',0,...
                      'epsilon_w',1,'tau',1,...
                      'epsilon_w_end',0.0,'wallAxis',2);
                  
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,...                      
                      'Dmu',-0.05);

    config = v2struct(optsPhys,optsNum);

    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium([],struct('solver','Newton'));
    %EX.ComputeEquilibrium(EX.optsPhys.rhoGas_sat);                     
    EX.IDC.plot(EX.GetRhoEq());
    EX.ComputeDynamics();
    res.fig_handles = EX.PlotDynamics();
end                 

