function Job_ComputeContactAngle_40()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[45,75],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',40);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34
 
    optsNum = struct('PhysArea',PhysArea,'FexNum',Fex_Num,...
                     'maxComp_y2',35,'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.375);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************        
    filename   = [dirData,filesep,'Job40_DisjoiningPressure_',getTimeStr(),'.txt'];       
    Struct2File(filename,config,['Computed at ',datestr(now)]);
        
  %  ConvergenceSurfaceTensions(config);
        
    close all;
    ChangeDirData();                
	
    CLT = ContactLine(config);     
    CLT.Preprocess();
    CLT.ComputeEquilibrium();  

    CLT.InitAnalysisGrid([-5 30],[0.5 18]);
    CLT.ComputeAdsorptionIsotherm('load'); %load \2DChebData\POF_FMT_ContactLine\deg90\IterativeContinuationPostProcess\2014_1_20_18_46
    CLT.PostProcess_2DDisjoiningPressure();
    
	%f1 = figure('Color','white','Position',[0 0 1000 600]);
    [f1,f2] = CLT.Post_HFrom2DDisjoiningPressure();
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces'],f1);
    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces.fig']);
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure'],f2);
    saveas(f2,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure.fig']);
    
    inset2(f1,f2,0.4,[0.26,0.55]);
    %inset2(f1,f2,0.35,[0.22,0.55]);
    close(f2);      
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure'],f1);
    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure.fig']);
        
    CLT.FittingAdsorptionIsotherm([10 14],1);
    CLT.SumRule_DisjoiningPotential();
    %***************************************  

end