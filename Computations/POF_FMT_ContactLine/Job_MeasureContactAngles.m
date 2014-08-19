function Job_MeasureContactAngles()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[35,65],...
                      'L1',5,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',90);
    
    V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',10,...
                     'V2Num',V2Num);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,...
                      'nSpecies',1,...
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
    
	%opts.epw_YCA = 0.7:0.005:1.6;
    opts.epw_YCA = 0.5:0.05:1.35;
    opts.config  = config;
    res1 = DataStorage('ContactAngleMeasurements',@MeasureYoungContactAngles,opts,[]);
    
    opts.epw_YCA = 1.35:0.005:1.51;
    opts.config  = config;
    res2 = DataStorage('ContactAngleMeasurements',@MeasureYoungContactAngles,opts,[]);
    
	f1 = figure('Color','white','Position',[0 0 800 800]);
    plot(res1.epw,180/pi*res1.theta_CA,'k','linewidth',1.5); hold on;    
    plot(res2.epw,180/pi*res2.theta_CA,'k','linewidth',1.5); hold on;    
    xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25); 
    set(gca,'fontsize',20); 
    xlim([0.5 1.52]);
    
    print2eps([dirData filesep 'ContactAngleMeasurements90'],gcf);
    saveas(gcf,[dirData filesep 'ContactAngleMeasurements90.fig']);   
    
    %opts.epw_YCA = 1.50:0.001:1.54;
    %opts.config  = config;
    %resG = DataStorage('ContactAngleMeasurements',@MeasureYoungContactAngles,opts,[]);
    
    %***********************************************************************
    %***********************************************************************
    %***********************************************************************
    
    close all;    

    opts90_a.config                            = config;    
    opts90_a.config.optsNum.PhysArea.alpha_deg = 90; 
    opts90_a.config.optsNum.maxComp_y2         = 15;
    opts90_a.config.optsNum.PhysArea.N         = [50,80];
    opts90_a.config.optsNum.PhysArea.L1        = 4; 
    opts90_a.config.optsNum.PhysArea.N2bound   = 14; 
    opts90_a.epw                               = 1.:0.02:1.08;
   % resM90_a = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts90_a,[]);    
    
    opts90_b     = opts90_a;
    opts90_b.epw = 1.1:0.02:1.16;
   % resM90_b = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts90_b,[]);    
    
    opts90_c     = opts90_a;
    opts90_c.epw = 0.55:0.05:1.;
    resM90_c = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts90_c,[]);
        
    opts90 = opts90_a;
    %opts90.config                            = config;
    opts90.config.optsNum.maxComp_y2         = 10; %TODO:15
    opts90.config.optsNum.PhysArea.alpha_deg = 90;   
    opts90.epw   = 1.2:0.02:1.36;
 %   resM90      = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts90,[]);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Redo computations on 60 [deg] grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	opts60.config                            = config;
    opts60.config.optsNum.PhysArea.alpha_deg = 60;
    opts60.config.optsNum.PhysArea.N         = [45,90];
    opts60.config.optsNum.PhysArea.L1        = 4;
    opts60.config.optsNum.PhysArea.N2bound   = 14;
    opts60.config.optsNum.maxComp_y2         = 15;
    opts60.epw                               = 1.2:0.02:1.3;%resM90.epw(abs(resM90.thetaM-40)<=10);
    
    %******************************
    
    resM60 = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts60,[]);
    
    opts60_b     = opts60;
    opts60_b.epw = 1.32:0.02:1.4;
    resM60_b = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts60_b,[]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Redo computations on 40 [deg] grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	opts40.config                            = config;
    opts40.config.optsNum.PhysArea.alpha_deg = 40;
    opts40.config.optsNum.PhysArea.N         = [45,80];
    opts40.config.optsNum.PhysArea.L1        = 4;     
    opts40.config.optsNum.PhysArea.N2bound   = 14;
    opts40.epw                               = 1.3:0.02:1.4;%resM90.epw(abs(resM90.thetaM-40)<=10);
    resM40 = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts40,[]);
    
    opts40_b                                   = opts40;
    opts40_b.config.optsNum.maxComp_y2         = 15;
    opts40_b.epw                               = 1.32:0.02:1.44;%resM90.epw(abs(resM90.thetaM-40)<=10);
    resM40_b = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts40_b,[]);
     
     
     
    
    close all;
    %*************************************************
    figure('Color','white','Position',[0 0 800 800]);
    plot(res.epw,180/pi*res.theta_CA,'k','linewidth',1.5); hold on;    
    %Filled: 15, unfilled:1
    %o - 90deg
    %s - 60deg (square)
    %d - 40deg
    Msize = 8;
    %**
%     %plot(resM90_a.epw,resM90_a.theta_YCA,'k-.');    
%     plot(resM90_a.epw,resM90_a.thetaM,'ko','MarkerFaceColor','k','MarkerSize',Msize);
%     %errorbar(resM90_a.epw,resM90_a.thetaM,resM90_a.thetaM_err,'k');
% 	%**
%     %plot(resM90_b.epw,resM90_b.theta_YCA,'k-.');    
%     plot(resM90_b.epw,resM90_b.thetaM,'ko','MarkerFaceColor','k','MarkerSize',Msize);
%     %errorbar(resM90_b.epw,resM90_b.thetaM,resM90_b.thetaM_err,'k');
%     %**
%     %plot(resM90.epw,resM90.theta_YCA,'k-.');
%     plot(resM90.epw,resM90.thetaM,'ko','MarkerSize',Msize);    
%     %errorbar(resM90.epw,resM90.thetaM,resM90.thetaM_err,'k');
%     %**
%      %plot(resM60.epw,resM60.theta_YCA,'g-.');
%      plot(resM60.epw,resM60.thetaM,'ks','MarkerFaceColor','k','MarkerSize',Msize);
%      plot(resM60_b.epw,resM60_b.thetaM,'ks','MarkerFaceColor','k','MarkerSize',Msize);
%      %errorbar(resM60.epw,resM60.thetaM,resM60.thetaM_err,'g');
% % 	%**
%      %plot(resM40.epw,resM40.theta_YCA,'m-.');
%      plot(resM40.epw,resM40.thetaM,'kd','MarkerSize',Msize);
%      %errorbar(resM40.epw,resM40.thetaM,resM40.thetaM_err,'m');    
%      
%      %plot(resM40_b.epw,resM40_b.theta_YCA,'m-.');
%      plot(resM40_b.epw,resM40_b.thetaM,'kd','MarkerFaceColor','k','MarkerSize',Msize);
%      %errorbar(resM40_b.epw,resM40_b.thetaM,resM40_b.thetaM_err,'m');         
%     
%     %plot(resM90.epw,resM90.thetaM,'ko');
%     %plot(resM40.epw,resM40.thetaM,'kx');
%     xlabel('$\varepsilon_w/\varepsilon$','Interpreter','Latex','fontsize',25);
%     ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',25);
%     set(gca,'fontsize',20);
%     %xlim([min(res.epw) max(res.epw)]);
    xlim([1 1.52]);
%     
    f1 = figure('Color','white','Position',[0 0 800 800]);
    plot(res.epw,180/pi*res.theta_CA,'k','linewidth',1.5); hold on;    
    xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25); 
    set(gca,'fontsize',20); 
    
    %************************************************
    f2 = figure('Color','white','Position',[0 0 800 800]);    
    thetaYCA = 180/pi*res.theta_CA;    
    
  %  plotErr(resM90,'o','w');
    plotErr(resM90_a,'o','k');
    plotErr(resM90_b,'o','k');
    plotErr(resM60,'s','k');
    plotErr(resM60_b,'s','k');
    %plotErr(resM40,'d','w');   
    plotErr(resM40_b,'d','k');   
    
    xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta -\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);
    
    inset2(f1,f2,0.45,[0.25,0.25]);
    close(f2);            
    
    str = [getTimeStr() , '_ContactAngles'];
    print2eps([dirData filesep 'ContactAngleMeasurements' filesep str],gcf);
    saveas(gcf,[dirData filesep 'ContactAngleMeasurements' filesep str '.fig']);
    
    
    %LDA: 
    %epsilonw - angle
    %1.5      - 88
    %2        - 60
    %2.4      - 33
    %2.57     - 15.5

    %config.optsNum.PhysArea.alpha_deg = degAngle;      
    
    
    function res = MeasureContactAngles(opts,h)
        
        epw     = opts.epw;
        
        CAM       = zeros(size(epw));
        err       = zeros(size(epw));
        theta_YCA = zeros(size(epw));
        
        configM = opts.config;               
        
        CL = ContactLineHS(configM);
        CL.Preprocess();
        
        for i = 1:length(epw)
            close all;
            %ChangeDirData();
            CL.SetV1(epw(i));
            CL.ComputeST();
            CL.ComputeEquilibrium();                        
            
          %  CL.PostProcess_FilmThickness([-10 0]);
            theta_YCA(i)    = CL.alpha_YCA*180/pi;
            [CAM(i),err(i)] = CL.MeasureContactAngle(2,[min(CL.optsNum.maxComp_y2-5) (CL.optsNum.maxComp_y2-1)]);
        end
        
        res.epw        = epw;
        res.thetaM     = CAM;
        res.thetaM_err = err;
        res.theta_YCA  = theta_YCA;
         
    end
    
    function res = MeasureYoungContactAngles(opts,h)

        epw_YCA  = opts.epw_YCA; 
        theta_CA = zeros(size(epw_YCA));
        omWG     = zeros(size(epw_YCA));
        omWL     = zeros(size(epw_YCA));
        
        %*********************************************************************
        %(1) Measure Contact Angles from Youngs Equation        
        opts.config.optsNum.maxComp_y2 = -1;        
        
        confP = opts.config;
        confP.optsNum.PhysArea.N2bound = 3;
        confP.optsNum.PhysArea.N = [60;3];
        CLP = ContactLineHS(confP);
        CLP.Preprocess();     
        omLG = CLP.ST_1D.om_LiqGas;            
        
        
        opts.config.optsNum.PhysArea.N = [1;100];
        CLN = ContactLineHS(opts.config);
        CLN.Preprocess();     
        
        for j = 1:length(epw_YCA)                  
            CLN.optsPhys.V1.epsilon_w = epw_YCA(j);      
            CLN.Compute1D('WL');            
            CLN.Compute1D('WG');
            
            omWG(j) = CLN.ST_1D.om_wallGas;
            omWL(j) = CLN.ST_1D.om_wallLiq;
            
            theta_CA(j) = ComputeContactAngle(omWG(j),omWL(j),omLG);
                        
        end 
        
        res.theta_CA = theta_CA;
        res.omWG     = omWG;
        res.omWL     = omWL;
        res.omLG     = omLG;
        res.epw      = epw_YCA;

        %*********************************************************************
        
    end

    function plotErr(resIn,str1,str2)
        mark = iseq(res.epw,resIn.epw);
        h = thetaYCA(mark);
        plot(res.epw(mark),(resIn.thetaM - h),[str1,'k'],'MarkerFaceColor',str2,'MarkerSize',Msize); hold on;
        %plot(h,(resIn.theta_YCA - h),['k',str2],'linewidth',1.5);  
    end

    function z = iseq(x,y)
        z = false(size(x));        
        for i = 1:length(x)
            if(sum(abs(y-x(i)) < 1e-6)>0)
                z(i) = true;
            end
        end
    end
end



%         PhysArea = struct('N',[40,40],'L1',2,'L2',2,'y2wall',0.,...
%                           'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);
% 
%         PhysArea.Conv  = struct('L',1,'L2',[],'N',[20,20]);
% 
%         Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
%                            'Ncircle',1,'N1disc',35,'N2disc',34);
% 
%         optsNum = struct('PhysArea',PhysArea,...
%                          'FexNum',Fex_Num,...
%                          'maxComp_y2',10,...
%                          'y1Shift',0,...c 
%                          'plotTimes_T',100,...
%                          'plotTimes_Interval',0.1);
% 
%         V1 = struct('V1DV1','Vext_Cart_7',...
%                             'epsilon_w',0.74,... %0.482218,...
%                             'epsilon_w_max',0.3,....
%                             'tau',20);
%                         
%         V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'R',1);
% 
%         optsPhys = struct('V1',V1,'V2',V2,...                      
%                           'kBT',0.75,...   
%                           'gamma',2,...
%                           'Inertial',1,...
%                           'Dmu',0.0,'nSpecies',1,...
%                           'sigmaS',1);      
%                       
%         configuration = v2struct(optsNum,optsPhys);
%         configName    = SaveConfig(configuration,'Configurations');