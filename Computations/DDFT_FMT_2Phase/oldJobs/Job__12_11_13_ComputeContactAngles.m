function Job__12_11_13_ComputeContactAngles()
%This Job should compute the contact lines for a range of contact angles, 
% such that a full study of e.g. the line tensions is made possible
% For enhanced accuracy, 60x60 points are used. 
% No dynamics is to be computed

    PhysArea = struct('N',[40,40],'L1_Skewed',2,'L2',2,'y2wall',0.,...
                      'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[20,20]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',35,'N2disc',34);

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',[],...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_Cart_7','epsilon_w',0.74);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'R',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);        
    kBT    = (0.75:0.02:0.85);
        
    config.optsNum.PhysArea.alpha_deg = 90;
    for i = 1:length(kBT)
        
        config.optsPhys.kBT  = kBT(i);
        
        config.optsNum.maxComp_y2         = 10;
        sol = ThreePhaseContactLine_FMT_BH(config,'FMT_ContactLine_Equilibrium_BH_40X40');                
        config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;                       
        
        config.optsNum.maxComp_y2         = 15;
        sol = ThreePhaseContactLine_FMT_BH(config,'FMT_ContactLine_Equilibrium_BH_40X40');        
        config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;                       
        
        config.optsNum.maxComp_y2         = 20;
        sol = ThreePhaseContactLine_FMT_BH(config,'FMT_ContactLine_Equilibrium_BH_40X40');
        config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;                       
        
        
    end
end