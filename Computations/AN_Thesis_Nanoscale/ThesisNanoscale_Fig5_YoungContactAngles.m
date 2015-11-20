function ThesisNanoscale_Fig5_YoungContactAngles()
    
    AddPaths('ThesisNanoscale');       
    
     config = ThesisNanoscale_GetStandardConfig(90);    
    
     %***************************************************
     %***************************************************
     % Get table of epw for certain contact angles
     %thY_V = [0.0,30,45,60,90,120,135,150,180];
     % epw_V = FindEpwFromContactAngle(config,thY_V);
     
     %****************************************************
     [res,f1] = ComputeYoungContactAngle(config,[0.45:0.05:1.0,1.0:0.005:1.24,1.24:0.001:1.29]); %[0.55:0.05:1.]     
     AddPaths('ThesisNanoscale');  
     
     pbaspect([1 1 1]); 
     ylim([0 135]);
     set(gcf,'Color','white','Position',[0 0 250 200]);
     fullName = SaveFigure('Fig1_YoungContactAngles');                             
    
    %****************************************************
    %****************************************************
    
    f2 = figure('Color','white','Position',[0 0 250 200]);    	
    thetaYCA = 180/pi*res.theta_CA;    
  
    %old
    %ComputeAndPlot(0.55:0.05:1.25,90,15,'o','k');    
    %ComputeAndPlot(0.45:0.05:0.8,120,15,'+','k');
    %ComputeAndPlot(1.3:0.05:1.44,40,15,'d','k');    
    %new    
    ComputeAndPlot(0.5:0.05:1.15,90,15,'o','k');    
    ComputeAndPlot(0.5:0.05:0.7,120,15,'+','k');
    ComputeAndPlot(1.0:0.05:1.15,60,15,'d','k');    
    %ComputeAndPlot(0.55:0.05:1.15,90,15,'o','k');    
  
    pbaspect([1 1 1]);
    xlim([0.5 1.2]);
    xlabel('$\LJWdepth$','Interpreter','Latex');
    ylabel('$\theta -\thYoung[^\circ]$','Interpreter','Latex');    
           
    %inset2(f1,f2,0.38,[0.3,0.25]); close(f2);  
    fullName = SaveFigure('Fig1_YoungContactAngles_Errors');
    
    config.optsPhys.V1.epsilon_w = 1.15;
    f3 = PlotContourLines(config);
    xlabel(''); ylabel('');
	f1 = hgload([fullName '.fig']);
    set(gcf,'Color','white','Position',[0 0 250 200]);    
    inset2(f1,f3,0.28,[0.32,0.22]);
    close(f3);
    
    ChangeDirData();    
    SaveFigure('Fig1_YoungContactAngles_Errors');
    
    function ComputeAndPlot(epw,alpha,maxY2,symbol,color)
        opts.config                            = config;                
        opts.config.optsNum.maxComp_y2         = maxY2; 
        opts.config.optsNum.PhysArea.alpha_deg = alpha;        
        opts.epw                               = epw;    
        
        resL = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts,[],[],{'config_optsPhys_V1_epsilon_w'});
        plotErr(resL,symbol,color);
    end
         
    function f3 = PlotContourLines(config)
        config.optsNum.PlotAreaCart = struct('y1Min',-15,'y1Max',5,...
                                             'y2Min',0.5,'y2Max',20,...%'zMax',4,...
                                             'N1',100,'N2',100);
        config.optsNum.maxComp_y2  = 15;                                    

        CL = ContactLineHS(config);
        CL.Preprocess();
        CL.ComputeEquilibrium(struct('solver','Picard'));     
        %CLT.ComputeEquilibrium();                                
        CL.PlotContourResults();
        
        y1L = [config.optsNum.PlotAreaCart.y1Min config.optsNum.PlotAreaCart.y1Max];
        close all;
        f3 = figure('Color','white','Position',[0 0 500 300]);                        
        %ha = area(y1L,[20,29],15,'FaceColor',0.8*[1 1 1]); hold on;
        CL.PlotContourResults({});
        plot(y1L,config.optsNum.maxComp_y2*[1 1],'k--','linewidth',1.5);        
        
        xlabel('$y_1$','Interpreter','Latex');
        ylabel('$y_2$','Interpreter','Latex');
        set(gca,'XTick',[-10 0]);
        
       % text('String','$y_{max}$','VerticalAlignment', 'top','HorizontalAlignment','center ','Position',[(y1L(1)+6) (config.optsNum.maxComp_y2 + 4)],'Interpreter','Latex');
    end    
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
            CL.ComputeEquilibrium(struct('solver','Picard'));     
            %CL.ComputeEquilibrium();                        
            
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
        %mark = iseq(res.epw,resIn.epw);
        %h    = thetaYCA(mark);
        for i = 1:length(resIn.epw)
            mark = iseq(res.epw,resIn.epw(i));            
            epw_h = res.epw(mark);
            th_h  = thetaYCA(mark);           
            if(sum(mark)==0)
                continue;
            end
            plot(epw_h(1),resIn.thetaM(i) - th_h(1),[str1,'k'],'MarkerFaceColor',str2,'MarkerSize',3); hold on;
        end
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