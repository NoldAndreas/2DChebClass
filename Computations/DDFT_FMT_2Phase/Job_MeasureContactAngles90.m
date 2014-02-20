function Job_MeasureContactAngles90()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw'],'ORG');    

    PhysArea = struct('N',[35,65],...
                      'L1_Skewed',5,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',90);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
%   Fex_Num   = struct('Fex','CarnahanStarling');
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',10,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************
    %Check convergence of surface tensions and errors of conact density
   % ConvergenceSurfaceTensions(config);
    
    %***********************************************************
    %Setup result file for this Job
    filename   = ([dirData filesep]); 
    filename   = [filename,'Job_MeasureContactAngles_epw_',...
                                getTimeStr(),'.txt'];
    
    %filename    = [dirData filesep subDir filesep 'Job__12_11_13_ComputeContactAngles_epw.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
    
    close all;    

    opts90_a.config                            = config;    
    opts90_a.config.optsNum.PhysArea.alpha_deg = 90; 
    opts90_a.config.optsNum.maxComp_y2         = 15;
    opts90_a.config.optsNum.PhysArea.N         = [50,80];
    opts90_a.config.optsNum.PhysArea.L1_Skewed = 4; 
    opts90_a.config.optsNum.PhysArea.N2bound   = 14; 
    opts90_a.epw                               = 1.:0.02:1.08;
    %resM90_a = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts90_a,[]);
    
    configT = opts90_a.config;
    configT.optsPhys.V1.epsilon_w = 1.37;
    configT.optsNum.maxComp_y2 = 6;
    CLT = ContactLine(configT);
    CLT.Preprocess();
    CLT.ComputeEquilibrium();    
    
end