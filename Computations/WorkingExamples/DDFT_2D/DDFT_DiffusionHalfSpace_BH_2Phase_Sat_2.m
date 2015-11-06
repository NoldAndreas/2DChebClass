function [config,res] = DDFT_DiffusionHalfSpace_BH_2Phase_Sat_2()
    
    Phys_Area = struct('shape','HalfSpaceSkewed','N',[30,30],...
                        'y2Min',0,'L1',4,'L2',10,'alpha',pi/2);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,'N2',100,...
                       'y2Min',0,'y2Max',20);
        
    Sub_Area = struct('shape','Box','N',[20,20],...
                      'y1Min',-1,'y1Max',1,...
                      'y2Min',0,'y2Max',2);    
                  	
    %V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[20,20]);    
    V2Num   = struct('Fex','SplitDisk','N',[20,20]);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:6,...
                     'V2Num',V2Num);                                           

    V1 = struct('V1DV1','Vext_Cart_7',...
                      'V0',0.0,'y10',2,'y20',0,...
                      'epsilon_w',1,'tau',1,'epsilon_w_end',0.0,...
                      'L1',4,'epsilon_w1',0.2,'w1_steepness',1); 
                  
    %V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,...                      
                      'Dmu',-0.03);

    config = v2struct(optsPhys,optsNum);

    AddPaths();
    EX     = DDFT_2D(config);
    
    EX.Preprocess(); 
    
%	yPtsCheck          = [20 10 ; 0 2 ; 0 0 ; -10 0 ];
%    EX.IDC.TestConvolutionMatrix(yPtsCheck,@Phi2DLongRange);
    EX.ComputeEquilibrium(EX.optsPhys.rhoGas_sat);                    
    EX.IDC.plot(EX.GetRhoEq());
    EX.ComputeDynamics();
    res.fig_handles = EX.PlotDynamics();
    
    function r = AnnMap(x,L)
        r = QuotientMap(x,L,1,inf);
    end
end                 

