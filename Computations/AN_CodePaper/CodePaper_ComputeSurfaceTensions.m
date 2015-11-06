function CodePaper_ComputeSurfaceTensions

    AddPaths('CodePaper');            
    close all;
    
    PhysArea = struct('N',[20,20],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                      
    V2Num    = struct('Fex','SplitDisk','N',[80,80]);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',FexNum,'V2Num',V2Num,...
                     'maxComp_y2',20,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.0);%1.375);%1.25)    

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);
                  
    config = v2struct(optsNum,optsPhys);          
    
%     epw    = [0.5:0.05:1.3];   
%     config.optsNum.V2Num.Fex     = 'ConstShortRange';
%     config.optsPhys.V2.V2DV2     = 'ConstShortRange'; 
%     config.optsPhys.V2.lambda    = 1.5;
%     config.optsPhys.V1.V1DV1     = 'Vext_SquareWell';    
%     [res{1},f1] = ComputeYoungContactAngle(config,epw);


    epw    = [0.5:0.05:1.3];   
  %  epw    = [0.85:0.01:0.86];   
    epw    = [1.15:0.01:1.16];   
	config.optsNum.V2Num.Fex     = 'SplitAnnulus';
	config.optsPhys.V2.V2DV2     = 'BarkerHenderson_2D';           
    config.optsPhys.V1.V1DV1     = 'Vext_BarkerHenderson_HardWall';    
    [res{1},f1] = ComputeYoungContactAngle(config,epw);
            
	epw    = [0.5:0.05:1.3];   
	config.optsNum.V2Num.Fex     = 'SplitAnnulus';
	config.optsPhys.V2.V2DV2     = 'BarkerHendersonCutoff_2D';           
    config.optsPhys.V1.V1DV1     = 'Vext_BarkerHenderson_HardWall';    
	[res{1},f1] = ComputeYoungContactAngle(config,epw);

    epw    = [1.2:0.05:1.6];   
    config.optsNum.V2Num.Fex     = 'SplitDisk';                
    config.optsPhys.V2.V2DV2     = 'ExponentialDouble';
    config.optsPhys.V1.V1DV1     = 'Vext_Exp_HardWall';
  	[res{2},f1] = ComputeYoungContactAngle(config,epw);
    
    epw    = [0.5:0.05:1.3];   
    config.optsNum.V2Num.Fex     = 'SplitDisk';
    config.optsPhys.V2.V2DV2     = 'BarkerHenderson_2D'; 
    config.optsPhys.V1.V1DV1     = 'Vext_BarkerHenderson_HardWall';    
	[res{3},f1] = ComputeYoungContactAngle(config,epw);                
    
end