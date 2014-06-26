function [optsNum,optsPhys] = Test_DDFT_EqBox_2Phase_2()
    
    Phys_Area = struct('shape','Box','N',[30,30],...
                       'y1Min',-4,'y1Max',4,...
                       'y2Min',-4,'y2Max',4);
    
    Plot_Area = Phys_Area;  Plot_Area.N1 = 100; Plot_Area.N2 = 100;
   
    FexNum       = struct('Fex','Meanfield');                          
    
    optsNum = struct('PhysArea',Phys_Area,'PlotArea',Plot_Area,...
                    'V2Num',FexNum);
                 
    V1 = struct('V1DV1','Vext_Cart_Capillary_Static',...                      
                      'V0',0.0,'y10',-5,'y20',0,'epsilon_w',[1,1,1,1]);
	V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,'Dmu',0);

    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium(EX.optsPhys.rhoGas_sat);                 

end                 