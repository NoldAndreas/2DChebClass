function Job_ComputeExactAdsorptionIsotherm()

    AddPaths();   
    global dirData    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[1,250],...
                      'L1',5,'L2',4,'L2_AD',4.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',90);

    V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.7);%1.375);%1.25);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.03,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    
    config.optsPhys.V1.epsilon_w = 0.7;%    1.0;%1.25;%0.55;% 1.375; %0.7;%1.25;%375;%25; %375;%47;%1.25;
    
    CL = ContactLineHS(config);
	CL.Preprocess();    
    CL.ComputeAdsorptionIsotherm(300,'drying');    
    %CL.ComputeAdsorptionIsotherm(300,'wetting');    
    CL.FittingAdsorptionIsotherm([10 14],1)
    if(optsPhys.kBT == 0.75)
        CL.SumRule_AdsorptionIsotherm(0.3463);
    end
end