function Job_ComputeContactAngle_20()  
    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[45,90],...
                      'L1_Skewed',4,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',20);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
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
    
    %try other iterative procedure 
    configT = config;
    configT.optsPhys.V1.epsilon_w = 1.47;
    configT.optsNum.maxComp_y2 = 15;
    
    %**************************************************
    %**************************************************
    %Plot Diagrams
    configT.optsNum.PhysArea.alpha_deg = 21;
    configT.optsNum.maxComp_y2 = 25;
    CLT = ContactLine(configT);     
    CLT.Preprocess();
    CLT.ComputeEquilibrium();  

    CLT.InitAnalysisGrid([-5 30],[0.5 18]);
    CLT.ComputeAdsorptionIsotherm('load'); %load 2014_1_30_16_37
    CLT.PostProcess_2DDisjoiningPressure();
    
    %f1 = figure('Color','white','Position',[0 0 1000 600]);
    [f1,f2] = CLT.Post_HFrom2DDisjoiningPressure();
    inset2(f1,f2,0.34,[0.28,0.5]);
    %inset2(f1,f2,0.35,[0.22,0.55]);
    % close(f2);      
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure'],f1);
    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure.fig']);
    
    CLT.FittingAdsorptionIsotherm([10 14],1);
    CLT.SumRule_DisjoiningPotential();
    %**************************************************
    %**************************************************
    
    CLT = ContactLine(configT);
    CLT.Preprocess();
    
    [y2_15,theta_15] = GetY2Theta(15);
    [y2_20,theta_20] = GetY2Theta(20);
    [y2_25,theta_25] = GetY2Theta(25);
    [y2_30,theta_30] = GetY2Theta(30);    
    
    configT.optsNum.PhysArea.alpha_deg = 21;
    configT.optsNum.maxComp_y2 = 15;
    CLT = ContactLine(configT);
    CLT.Preprocess();
    
    
    [y2_21_15,theta_21_15] = GetY2Theta(15);
    [y2_21_20,theta_21_20] = GetY2Theta(20);
    [y2_21_25,theta_21_25] = GetY2Theta(25);
%    [y2_21_30,theta_21_30] = GetY2Theta(30);    
        
    f1 = figure('Color','white','Position',[0 0 800 800]);    
    plot(y2_15,theta_15,'k:','linewidth',1.5); hold on; 
    plot(y2_20,theta_20,'k-.','linewidth',1.5); hold on; 
    plot(y2_25,theta_25,'k--','linewidth',1.5); hold on; 
    plot(y2_30,theta_30,'k','linewidth',1.5); hold on; 
    
    plot(y2_21_15,theta_21_15,'b:','linewidth',1.5); hold on; 
    plot(y2_21_20,theta_21_20,'b-.','linewidth',1.5); hold on; 
    plot(y2_21_25,theta_21_25,'b--','linewidth',1.5); hold on; 
    %plot(y2_21_30,theta_21_30,'b','linewidth',1.5); hold on; 
    
    xlim([7 22]);
    ylim([20.5 24.5]);
    
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'linewidth',1.5,'fontsize',20);    
    
    print2eps([dirData filesep 'CA_Asymptotics' filesep 'CA20'],f1);
    saveas(f1,[dirData filesep 'CA_Asymptotics' filesep 'CA20.fig']);

    
   function [y2,theta] = GetY2Theta(y2Max)
        CLT.optsNum.maxComp_y2 = y2Max;
        CLT.ComputeEquilibrium();  
        [y2,theta] = CLT.PlotInterfaceAnalysisY2([5 (y2Max+3)]);
   end 
  
end