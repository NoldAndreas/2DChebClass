function MMNP_Fig1_YoungContactAngles()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'MMNP'],'ORG');    

    PhysArea = struct('N',[50,80],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,'N2bound',14,'h',1,...
                      'alpha_deg',90);
     
    V2Num     = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
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
        
    [res,f1] = ComputeYoungContactAngle(config,[0.45:0.05:1.35,1.35:0.005:1.51]);
    %[res,f1] = ComputeYoungContactAngle(config,[0.55:0.05:1.]);
    
    SaveFigure(['Figures' filesep 'Fig1']);
            
    f2 = figure('Color','white','Position',[0 0 800 800]);    	
    thetaYCA = 180/pi*res.theta_CA;    
    
    %ComputeAndPlot(0.55:0.02:1.25,90,15,'o','k');    
    ComputeAndPlot(0.55:0.05:1.25,90,15,'o','k');    
    ComputeAndPlot(0.45:0.05:0.8,120,15,'+','k');
    ComputeAndPlot(1.3:0.05:1.44,40,15,'d','k');    
  %  ComputeAndPlot(1.2:0.05:1.4,60,15,'s','k'); % N = [45,90]?? 
    
    xlim([0.45 1.44]);
    xlabel('${\alpha_w \sigma^3}/{\varepsilon}$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta -\theta_{Y}[^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);    
           
    inset2(f1,f2,0.38,[0.3,0.25]);
    close(f2);  
    fullName = SaveFigure('temp');
   
    config.optsPhys.V1.epsilon_w = 0.55;
    f3 = PlotContourLines(config);
	f1 = hgload([fullName '.fig']);
    inset2(f1,f3,0.35,[0.55,0.65]);
    close(f3);
    
    ChangeDirData();    
    SaveFigure(['Figures' filesep 'Fig1']);
         
    function f3 = PlotContourLines(config)
        config.optsNum.PlotAreaCart = struct('y1Min',-4,'y1Max',20,...
                                             'y2Min',0.5,'y2Max',20,...%'zMax',4,...
                                             'N1',100,'N2',100);
        
        CL = ContactLineHS(config);
        CL.Preprocess();
        CL.ComputeEquilibrium();                                
        CL.PlotContourResults();
        
        close all;
        f3 = figure('Color','white','Position',[0 0 500 300]);                        
        ha = area([-7,20],[20,29],15,'FaceColor',0.8*[1 1 1]); hold on;
        CL.PlotContourResults({});
        plot([-5 20],[15 15],'k--','linewidth',2);
        set(gca,'linewidth',2);        
        text('String','$y_{max}$','VerticalAlignment', 'top','HorizontalAlignment','center ','Position',[12 17],'Interpreter','Latex','fontsize',20);        
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
        %mark = iseq(res.epw,resIn.epw);
        %h    = thetaYCA(mark);
        for i = 1:length(resIn.epw)
            mark = iseq(res.epw,resIn.epw(i));
            epw_h = res.epw(mark);
            th_h  = thetaYCA(mark);
            plot(epw_h(1),resIn.thetaM(i) - th_h(1),[str1,'k'],'MarkerFaceColor',str2,'MarkerSize',8); hold on;
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
    function ComputeAndPlot(epw,alpha,maxY2,symbol,color)
        opts.config = config;                
        opts.config.optsNum.maxComp_y2         = maxY2; 
        opts.config.optsNum.PhysArea.alpha_deg = alpha;        
        opts.epw    = epw;    
        
        resL = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts,[]);
        plotErr(resL,symbol,color);
    end
end