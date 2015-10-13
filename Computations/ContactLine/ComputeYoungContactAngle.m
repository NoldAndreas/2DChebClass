function [res,f1] = ComputeYoungContactAngle(config,epw)
    %Check convergence of surface tensions and errors of conact density
   % ConvergenceSurfaceTensions(config);
    
    %***********************************************************
    %Setup result file for this Job
    global dirData
    
    filename   = ([dirData filesep]); 
    filename   = [filename,'Job_MeasureContactAngles_epw_',...
                                getTimeStr(),'.txt'];
    
    %filename    = [dirData filesep subDir filesep 'Job__12_11_13_ComputeContactAngles_epw.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
    	
    opts.epw_YCA = epw;
    opts.config  = config;
    res          = DataStorage('YoungContactAngleMeasurements',@MeasureYoungContactAngles,opts,[],[],{'config_optsNum_PhysArea_N'});
    
  	f1 = figure('Color','white','Position',[0 0 300 250]);
    plot(res.epw,180/pi*res.theta_CA,'k','linewidth',1.0); hold on;        
    xlabel('$\LJWdepth$','Interpreter','Latex');
    ylabel('$\thYoung[^\circ]$','Interpreter','Latex'); 
    set(gca,'fontsize',20); 
    xlim([0.5 1.3]);
    set(gca,'YTick',[0 45 90 135]);
    
    SaveFigure('ContactAngleMeasurements90');
    
   
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
        confP.optsNum.PhysArea.N       = [60;3];
        CLP = ContactLineHS(confP);
        CLP.Preprocess();     close all;
        omLG = CLP.ST_1D.om_LiqGas;            
        
        
        opts.config.optsNum.PhysArea.N = [1;100];
        CLN = ContactLineHS(opts.config);
        CLN.Preprocess();     close all;
        
        for j = 1:length(epw_YCA)                  
            CLN.optsPhys.V1.epsilon_w = epw_YCA(j);      
            CLN.Compute1D('WL'); close all;
            CLN.Compute1D('WG'); close all;
            
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
end