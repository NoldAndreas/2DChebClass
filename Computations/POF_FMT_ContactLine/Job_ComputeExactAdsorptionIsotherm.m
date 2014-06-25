function Job_ComputeExactAdsorptionIsotherm()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    
    %FMT_CLEq_BH_40X40_epw

    PhysArea = struct('N',[1,250],...
                      'L1',5,'L2',4,'L2_AD',4.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',90);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
%   Fex_Num   = struct('Fex','CarnahanStarling');
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.375);%1.25);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************
    %Check convergence of surface tensions and errors of conact density
    
    config.optsPhys.V1.epsilon_w = 1.25; %375;%47;%1.25;
    
    CL = ContactLine(config); %CL with normal high resolution
	CL.Preprocess();     
    CL.ComputeAdsorptionIsotherm(300);
    
    CL.FittingAdsorptionIsotherm([10 14],1)
end