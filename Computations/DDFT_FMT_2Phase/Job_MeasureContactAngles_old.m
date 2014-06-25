function Job_MeasureContactAngles()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw'],'ORG');    

    PhysArea = struct('N',[35,65],...
                      'L1',5,'L2',2,'L2_AD',2.,...
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
    %Setup result file for this Job
    filename   = ([dirData filesep]); 
    filename   = [filename,'Job_MeasureContactAngles_epw_',...
                                getTimeStr(),'.txt'];
    
    %filename    = [dirData filesep subDir filesep 'Job__12_11_13_ComputeContactAngles_epw.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
    
%	opts.epw_YCA = 0.7:0.005:1.6;
%    opts.config  = config;
%    res = DataStorage('ContactAngleMeasurements',@MeasureYoungContactAngles,opts,[]);
%    close all;
    
 %   optss.config                    = config;
 %   optss.config.optsNum.maxComp_y2 = 10;
 %   optss.config.optsNum.PhysArea.alpha_deg = 90;
 %   optss.epw                      = 1.2:0.02:1.4;
 %   resM90 = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,optss,[]);
      
    close all;
      
	opts40.config                            = config;
    opts40.config.optsNum.PhysArea.alpha_deg = 40;
    opts40.config.optsNum.PhysArea.N         = [45,80];
    opts40.config.optsNum.PhysArea.L1        = 4;     
    opts40.config.optsNum.PhysArea.N2bound   = 14;
    opts40.epw                               = 1.3:0.02:1.4;%sresM90.epw(abs(resM90.thetaM-40)<=10);
    resM40 = DataStorage('ContactAngleMeasurements',@MeasureContactAngles,opts40,[]);
    
    figure('Color','white','Position',[0 0 800 800]);
    plot(res.epw,180/pi*res.theta_CA,'k','linewidth',1.5); hold on;
    plot(resM90.epw,resM90.thetaM,'ko');
    plot(resM40.epw,resM40.thetaM,'kx');
    xlabel('$\varepsilon_w/\varepsilon$','Interpreter','Latex','fontsize',20);
    ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);
    xlim([min(res.epw) max(res.epw)]);
            
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
        CAM     = zeros(size(epw));
        configM = opts.config;               
        
        CL = ContactLine(configM);
        CL.Preprocess();
        
        for i = 1:length(epw)
            close all;
            %ChangeDirData();

            CL.optsPhys.V1.epsilon_w  = epw(i); 
            CL.ComputeST();
            CL.ComputeEquilibrium();
            CL.PostProcess_FilmThickness([-10 0]);
            CAM(i) = CL.CA_deg_measured;
        end
        
        res.epw    = epw;
        res.thetaM = CAM;
         
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
        CLP = ContactLine(confP);
        CLP.Preprocess();     
        omLG = CLP.ST_1D.om_LiqGas;            
        
        
        opts.config.optsNum.PhysArea.N = [1;100];
        CLN = ContactLine(opts.config);
        CLN.Preprocess();     
        
        for j = 1:length(epw_YCA)            
            CLN.optsPhys.V1.epsilon_w = epw_YCA(j);                        
            CLN.ComputeST(false);
            
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