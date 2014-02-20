function Job_ComputeContactAngle_66()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw'],'ORG');    

    PhysArea = struct('N',[35,65],...
                      'L1_Skewed',4,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',66);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34
 
    optsNum = struct('PhysArea',PhysArea,'FexNum',Fex_Num,...
                     'maxComp_y2',35,'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.2);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************        
    filename   = [dirData,filesep,'Job66_DisjoiningPressure_',getTimeStr(),'.txt'];       
    Struct2File(filename,config,['Computed at ',datestr(now)]);
        
    ConvergenceSurfaceTensions(config);
        
    close all;
    ChangeDirData();                
	CL = ContactLine(config);
	CL.Compute();

end