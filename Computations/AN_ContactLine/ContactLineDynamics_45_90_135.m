function ContactLineDynamics_45_90_135(kBT)
    
    if(nargin < 1)
        kBT = 0.75;
    end

    PhysArea = struct('N',[40,40],...
                      'L1',4,'L2',2,...                        
                      'alpha_deg',[]);
                  
	SubArea      = struct('shape','Box','y1Min',-2,'y1Max',2,...
                          'y2Min',0.5,'y2Max',2.5,...
                          'N',[40,40]);
                              
                      
    V2Num    = struct('Fex','SplitAnnulus','N',[60,60]);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',40,'N2disc',40);
                   
	plotTimes = struct('t_int',[0,400],'t_n',100);

    optsNum = struct('PhysArea',PhysArea,... %'PlotAreaCart',PlotAreaCart,...                     
                     'FexNum',FexNum,'V2Num',V2Num,...
                     'SubArea',SubArea,...
                     'maxComp_y2',10,...%'y1Shift',0,...
                     'plotTimes',plotTimes);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',[]);%,...%1.154 = 45 degrees
                %'tau',5,'epsilon_w_end',1.0);
            
    %optsViscosity = struct('etaC',1,'zetaC',0);    
    %optsViscosity = struct('etaL1',2,'zetaC',1);
    
    
    %******************************************************
    %******************************************************
    %Test on 11/03/2016
	optsViscosity = struct('etaLiq',3,'etaVap',0.1,...
                           'zetaLiq',3,'zetaVap',0.1);                       
    %previous version
    %optsViscosity = struct('etaLiq',5,'etaVap',1,...
    %                           'zetaLiq',5,'zetaVap',1);
    %******************************************************
    %******************************************************
                       
    %BCwall        = struct('bc','sinHalf','tau',1);
	%BCwall        = struct('bc','exp','tau',1,'u_max',0.2);
    BCwall = [];

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1,...
                      'Inertial',true,'gammaS',0,...%4% 'Fext',[-1,0],...%);%,...% %'BCWall_U',BCwall,...%);                     
                      'viscosity',optsViscosity);	

    config = v2struct(optsNum,optsPhys);      

    %config.optsPhys.kBT = 0.9;            
    if(config.optsPhys.kBT == 0.9)
        config.optsNum.plotTimes.t_int(2) = 800;
        config.optsNum.PhysArea.L1 = 6;
    end
     
    ContactLineDynamics_X_degrees(config,{'90','advancing','snapshots'});
    %ContactLineDynamics_X_degrees(config,{'45','advancing','snapshots'});        
    %ContactLineDynamics_X_degrees(config,{'135','advancing','snapshots'});
    
    config.optsNum.plotTimes.t_int(2) = 200;    
    %ContactLineDynamics_X_degrees(config,{'90','receding','snapshots'});   
    ContactLineDynamics_X_degrees(config,{'45','receding','snapshots'});        
    ContactLineDynamics_X_degrees(config,{'135','receding','snapshots'});
    
end