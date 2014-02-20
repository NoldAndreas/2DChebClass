function Job_MeasureContactAngles60()  
    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

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
    
    %opts.epw_YCA = 1.50:0.001:1.54;
    %opts.config  = config;
    %resG = DataStorage('ContactAngleMeasurements',@MeasureYoungContactAngles,opts,[]);
   
    close all;    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Redo computations on 60 [deg] grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	opts60.config                            = config;
    opts60.config.optsNum.PhysArea.alpha_deg = 60;
    opts60.config.optsNum.PhysArea.N         = [45,90];
    opts60.config.optsNum.PhysArea.L1_Skewed = 4;
    opts60.config.optsNum.PhysArea.N2bound   = 14;
    opts60.config.optsNum.maxComp_y2         = 15;
    opts60.epw                               = 1.2:0.02:1.3;%resM90.epw(abs(resM90.thetaM-40)<=10);
    
    %******************************
    %try other iterative procedure
    configT = opts60.config;
    configT.optsPhys.V1.epsilon_w = 1.25;
    configT.optsNum.maxComp_y2 = 20;
    CLT = ContactLine(configT);
    CLT.Preprocess();
    %CLT.ComputeEquilibrium();    
    %**************************************        
    
    [y2_20,theta_20] = GetY2Theta(20);
    [y2_25,theta_25] = GetY2Theta(25);
    [y2_30,theta_30] = GetY2Theta(30);
    [y2_35,theta_35] = GetY2Theta(35);
    
    f1 = figure('Color','white','Position',[0 0 800 800]);    
    plot(y2_20,theta_20,'k:','linewidth',1.5); hold on; 
    plot(y2_25,theta_25,'k-.','linewidth',1.5); hold on; 
    plot(y2_30,theta_30,'k--','linewidth',1.5); hold on; 
    plot(y2_35,theta_35,'k','linewidth',1.5); hold on; 
    
    xlim([10 37]);
    ylim([60 60.8]);   
    set(gca,'YTick',60:0.2:60.8);
    
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'linewidth',1.5,'fontsize',20);    

    f2 = figure('Color','white');
    CLT.PlotEquilibriumResults([-2 12],[0.5 14.5],true,false);   
        
    inset2(f1,f2,0.7,[0.35,0.45]);
    close(f2);   
    
    print2eps([dirData filesep 'CA_Asymptotics' filesep 'CA60'],f1);
    saveas(f1,[dirData filesep 'CA_Asymptotics' filesep 'CA60.fig']);
        
    %**************************************
     
%     configT.optsNum.maxComp_y2 = 25;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();  
%     
%     configT.optsNum.maxComp_y2 = 30;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium(); 
%     
%     configT.optsNum.maxComp_y2 = 35;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium(); 
    
	function [y2,theta] = GetY2Theta(y2Max)
        CLT.optsNum.maxComp_y2 = y2Max;
        CLT.ComputeEquilibrium();  
        [y2,theta] = CLT.PlotInterfaceAnalysisY2([5 (y2Max+3)]);        
    end
end