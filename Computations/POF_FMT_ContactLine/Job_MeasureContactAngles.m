function Job_MeasureContactAngles()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[50,80],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,'N2bound',14,'h',1,...
                      'alpha_deg',90);
    
    V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',15,...
                     'V2Num',V2Num);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,...
                      'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
        
    [res,f1] = ComputeYoungContactAngle(config,[0.5:0.05:1.35,1.35:0.005:1.51]);
    %[res,f1] = ComputeYoungContactAngle(config,[0.55:0.05:1.]);
            
    f2 = figure('Color','white','Position',[0 0 800 800]);    	
    thetaYCA = 180/pi*res.theta_CA;    
    
    ComputeAndPlot(0.55:0.02:1.25,90,15,'o','k');    
    %ComputeAndPlot(0.55:0.05:1.,90,15,'o','k');    
    ComputeAndPlot(1.2:0.02:1.4,60,15,'s','k'); % N = [45,90]??        
    ComputeAndPlot(1.3:0.02:1.44,40,15,'o','k');    
    
    xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta -\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);    
           
    inset2(f1,f2,0.45,[0.5,0.5]);
    close(f2);            
    
    str = [getTimeStr() , '_ContactAngles'];
    print2eps([dirData filesep 'ContactAngleMeasurements' filesep str],gcf);
    saveas(gcf,[dirData filesep 'ContactAngleMeasurements' filesep str '.fig']);
    
    close all;        
    %*************************************************
%     figure('Color','white','Position',[0 0 800 800]);
% 	f1 = figure('Color','white','Position',[0 0 800 800]);
%     plot(res.epw,180/pi*res.theta_CA,'k','linewidth',1.5); hold on;    
%     xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
%     ylabel('$\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25); 
%     set(gca,'fontsize',20);
        
    %Filled: 15, unfilled:1
    %o - 90deg
    %s - 60deg (square)
    %d - 40deg
   % Msize = 8;
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
%    xlim([1 1.52]);
%        

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

    function plotErr(resIn,str1,str2)
        mark = iseq(res.epw,resIn.epw);
        h = thetaYCA(mark);
        plot(res.epw(mark),(resIn.thetaM - h),[str1,'k'],'MarkerFaceColor',str2,'MarkerSize',8); hold on;
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

  function ComputeAndPlot(epw,alpha,maxY2,symbol,color)
        opts.config = config;                
        opts.config.optsNum.maxComp_y2         = maxY2; 
        opts.config.optsNum.PhysArea.alpha_deg = alpha;        
        opts.epw    = epw;    
        
        resL = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts,[]);
        plotErr(resL,symbol,color);
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