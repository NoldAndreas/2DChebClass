function Job_MeasureContactAngles40()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw'],'ORG');    

    PhysArea = struct('N',[35,65],...
                      'L1',5,...
                      'L2',2,...
                      'L2_AD',2.,...
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
    %Redo computations on 40 [deg] grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	opts40.config                            = config;
    opts40.config.optsNum.PhysArea.alpha_deg = 42.25 ;
    opts40.config.optsNum.PhysArea.N         = [45,80];
    opts40.config.optsNum.PhysArea.L1        = 4;     
    opts40.config.optsNum.PhysArea.N2bound   = 14;
    opts40.epw                               = 1.3:0.02:1.4;%resM90.epw(abs(resM90.thetaM-40)<=10);
    %resM40 = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts40,[]);
    
    
    configT = opts40.config;
    configT.optsPhys.V1.epsilon_w = 1.37;    
    configT.optsNum.PhysArea.alpha_deg = 40;%42.25 ;
    CLT = ContactLine(configT);
    %CLT.Preprocess();    
    res = DataStorage('CA_Asymptotics',@ComputeAll_Asymptotics,configT,[]);
    
    close all;
    f1 = figure('Color','white','Position',[0 0 1000 800]);
    plotGraphs(res.CN_40,'r',1.5);
    plotGraphs(res.CN_42,'b',1.5);  
    plotGraphs(res.CN_4225,'k',1.5);
    %plotGraphs(res.CN_425,'g',1.5);
    xlim([7 17]);
    ylim([42 45]);
    
    xlabel('$y$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta$','Interpreter','Latex','fontsize',25);
    set(gca,'linewidth',1.5,'fontsize',20);    
    
    print2eps([dirData filesep 'CA_Asymptotics' filesep 'CA40'],f1);
    saveas(f1,[dirData filesep 'CA_Asymptotics' filesep 'CA40.fig']);

    
    function res = ComputeAll_Asymptotics(configT,h)
               
        configT.optsNum.PhysArea.alpha_deg = 40;%42.25 ;
        CLT = ContactLine(configT);
        CLT.Preprocess();

        [res.CN_40.y2(:,1),res.CN_40.theta(:,1)] = GetY2Theta(15);
        [res.CN_40.y2(:,2),res.CN_40.theta(:,2)] = GetY2Theta(20);
        [res.CN_40.y2(:,3),res.CN_40.theta(:,3)] = GetY2Theta(25);
        [res.CN_40.y2(:,4),res.CN_40.theta(:,4)] = GetY2Theta(30);                
        [res.CN_40.y2(:,5),res.CN_40.theta(:,5)] = GetY2Theta(35);

        configT.optsNum.PhysArea.alpha_deg = 42;%42.25 ;
        CLT = ContactLine(configT);
        CLT.Preprocess();
        
        [res.CN_42.y2(:,1),res.CN_42.theta(:,1)] = GetY2Theta(15);
        [res.CN_42.y2(:,2),res.CN_42.theta(:,2)] = GetY2Theta(20);
        [res.CN_42.y2(:,3),res.CN_42.theta(:,3)] = GetY2Theta(25);
        [res.CN_42.y2(:,4),res.CN_42.theta(:,4)] = GetY2Theta(30);
        [res.CN_42.y2(:,5),res.CN_42.theta(:,5)] = GetY2Theta(35);

        configT.optsNum.PhysArea.alpha_deg = 42.5;%42.25 ;
        CLT = ContactLine(configT);
        CLT.Preprocess();
        
        [res.CN_425.y2(:,1),res.CN_425.theta(:,1)] = GetY2Theta(15);
        [res.CN_425.y2(:,2),res.CN_425.theta(:,2)] = GetY2Theta(20);
        [res.CN_425.y2(:,3),res.CN_425.theta(:,3)] = GetY2Theta(25);
        [res.CN_425.y2(:,4),res.CN_425.theta(:,4)] = GetY2Theta(30);                
        [res.CN_425.y2(:,5),res.CN_425.theta(:,5)] = GetY2Theta(35);
        
        configT.optsNum.PhysArea.alpha_deg = 42.25;%42.25 ;
        CLT = ContactLine(configT);
        CLT.Preprocess();

        [res.CN_4225.y2(:,1),res.CN_4225.theta(:,1)] = GetY2Theta(15);
        [res.CN_4225.y2(:,2),res.CN_4225.theta(:,2)] = GetY2Theta(20);
        [res.CN_4225.y2(:,3),res.CN_4225.theta(:,3)] = GetY2Theta(25);
        [res.CN_4225.y2(:,4),res.CN_4225.theta(:,4)] = GetY2Theta(30);    
        [res.CN_4225.y2(:,5),res.CN_4225.theta(:,5)] = GetY2Theta(35);    
        
    end
    
    function [y2,theta] = GetY2Theta(y2Max)
        CLT.optsNum.maxComp_y2 = y2Max;
        CLT.ComputeEquilibrium();  
        [y2,theta] = CLT.PlotInterfaceAnalysisY2([5 (y2Max+3)]);
    end

    function plotGraphs(res,c,lw)        
        plot(res.y2(:,1),res.theta(:,1),[c,':'],'linewidth',lw); hold on;
        plot(res.y2(:,2),res.theta(:,2),[c,'-.'],'linewidth',lw); hold on;
        plot(res.y2(:,3),res.theta(:,3),[c,'--'],'linewidth',lw); hold on;
        plot(res.y2(:,4),res.theta(:,4),[c,'-'],'linewidth',lw); hold on;
       % plot(res.y2(:,5),res.theta(:,5),[c,'-']); hold on;
    end
    
end

%     
%     CLT.optsNum.maxComp_y2 = 20;
%     CLT.ComputeEquilibrium();  
%     [y2_40_15,theta_40_15] = CLT.PlotInterfaceAnalysisY2([5 18]);
%     
%     CLT.optsNum.maxComp_y2 = 15;
%     CLT.ComputeEquilibrium();  
%     [y2_40_15,theta_40_15] = CLT.PlotInterfaceAnalysisY2([5 18]);
% 
%     configT = opts40.config;
%     configT.optsPhys.V1.epsilon_w = 1.37;
%     configT.optsNum.maxComp_y2 = 20;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();    
%     [y2_40_20,theta_40_20] = CLT.PlotInterfaceAnalysisY2([5 23]);
% 
%     configT = opts40.config;
%     configT.optsPhys.V1.epsilon_w = 1.37;
%     configT.optsNum.maxComp_y2 = 25;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();    
%     [y2_40_20,theta_40_20] = CLT.PlotInterfaceAnalysisY2([5 23]);    
% 
%     configT = opts40.config;
%     configT.optsPhys.V1.epsilon_w = 1.37;
%     configT.optsNum.maxComp_y2 = 30;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();  
%     CLT.PlotInterfaceAnalysisY2([10 38]);
%     
%     configT = opts40.config;
%     configT.optsPhys.V1.epsilon_w = 1.37;
%     configT.optsNum.maxComp_y2 = 35;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();    
%     CLT.PlotInterfaceAnalysisY2([10 38]);

     