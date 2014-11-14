function [optsNum,optsPhys] = DFT_EqInfCapillary_Slit()

    %Numerical Parameters    
    Phys_Area = struct('shape','InfCapillary','N',[25,20],...
                       'y1Min',-inf,'y1Max',inf,'L1',5,...
                       'y2Min',-4,'y2Max',4);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'L1',5,'N1',100,...
                       'y2Min',-4,'y2Max',4,'N2',100);
                   
    FexNum       = struct('Fex','Meanfield','L1',3,'L2',3,'N',[30,30]);                          
            
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',FexNum);
                     %'DDFTCode','DDFT_EqBoxInf_1Phase','Cutoff',5);                            
                 
    V1 = struct('V1DV1','Vext_Cart_Slit_Static','epsilon_w',[1,1,1,1]);
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,'nParticlesS',0);%'Dmu',0);
                      %'V2DV2_1D','Phi1DLongRange');                 

    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    
    rL = EX.optsPhys.rhoLiq_sat;
    rG = EX.optsPhys.rhoGas_sat;
    rhoIG = (rL+rG)/2 + (rL-rG)/2*tanh(EX.IDC.GetCartPts.y1_kv);
    
    EX.ComputeEquilibrium(rhoIG);
    EX.IDC.doPlots(EX.GetRhoEq());
end                 

