function Job_ComputeContactAngle_60()

    global dirData
    AddPaths();
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    %'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[45,90],...
                      'L1_Skewed',4,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',60);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34
 
    optsNum = struct('PhysArea',PhysArea,'FexNum',Fex_Num,...
                     'maxComp_y2',30,'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.25);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************        
    filename   = [dirData,filesep,'Job60_DisjoiningPressure_',getTimeStr(),'.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
        
  %  ConvergenceSurfaceTensions(config);
   
    
    %******************************
    %try other iterative procedure
   
    close all;
    ChangeDirData();                
    CLT = ContactLine(config);     
    CLT.Preprocess();
    CLT.ComputeEquilibrium();  

    CLT.optsNum.PlotArea.y1Min = -5;
    CLT.optsNum.PlotArea.y1Max = 30;
    
    CLT.InitAnalysisGrid([-5 20],[0.5 30]);%[-5 15],[0.5 18]);
    CLT.ComputeAdsorptionIsotherm(); %load 2014_1_22_10_9 
    CLT.PostProcess_2DDisjoiningPressure();

    [f1,f2] = CLT.Post_HFrom2DDisjoiningPressure();
    %inset2(f1,f2,0.4,[0.26,0.55]);
    set(gca,'XAxisLocation','top');
    inset2(f1,f2,0.35,[0.62,0.28]);
    close(f2);      
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure'],f1);
    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure.fig']);
        
        
    CLT.FittingAdsorptionIsotherm([10 14],1);
    CLT.SumRule_DisjoiningPotential();

end