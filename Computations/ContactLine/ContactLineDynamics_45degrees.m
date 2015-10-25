function ContactLineDynamics_45degrees()

    AddPaths('CodePaper');            
    close all;
    
    advancing = false;
    if(~advancing)
        alpha_deg = 45;
        epw       = 0.856; %= 90 degree contact angle
        maxT      = 400;
    else
        alpha_deg = 60;
        epw       = 1.22;  %= +/- 30 degree contact angle
        %epw       = 1.154; %= 45 degree contact angle
        maxT      = 400;
    end
    
    PhysArea = struct('N',[40,40],...
                      'L1',4,'L2',2,...                        
                      'alpha_deg',alpha_deg);
                  
	SubArea      = struct('shape','Box','y1Min',-2,'y1Max',2,...
                          'y2Min',0.5,'y2Max',2.5,...
                          'N',[40,40]);
                          
%     PlotAreaCart =     struct('y1Min',-5,'y1Max',20,...
%                               'y2Min',0.5,'y2Max',15.5,...
%                               'N1',100,'N2',100,'NFlux',40);

     PlotAreaCart =     struct('y1Min',-5,'y1Max',10,...
                               'y2Min',0.5,'y2Max',15.5,...
                               'N1',100,'N2',100,'NFlux',20);
                      
    V2Num    = struct('Fex','SplitAnnulus','N',[80,80]);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);
                   
	plotTimes = struct('t_int',[0,maxT],'t_n',100);

    optsNum = struct('PhysArea',PhysArea,...
                     'PlotAreaCart',PlotAreaCart,...
                     'FexNum',FexNum,'V2Num',V2Num,...
                     'SubArea',SubArea,...
                     'maxComp_y2',10,...%'y1Shift',0,...
                     'plotTimes',plotTimes);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);%,...%1.154 = 45 degrees
                %'tau',5,'epsilon_w_end',1.0);
            
    %optsViscosity = struct('etaC',1,'zetaC',0);    
    %optsViscosity = struct('etaL1',2,'zetaC',1);
    optsViscosity = struct('etaLiq',5,'etaVap',1,...
                           'zetaLiq',5,'zetaVap',1);
                            %'zetaC',1);
    %BCwall        = struct('bc','sinHalf','tau',1);
	BCwall        = struct('bc','exp','tau',1,'u_max',0.2);

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1,...
                      'Inertial',true,'gammaS',0,...%4% 'Fext',[-1,0],...%);%,...% 'BCWall_U',BCwall,...%);                     
                      'viscosity',optsViscosity);	

    config = v2struct(optsNum,optsPhys);                                    
        
    CL = ContactLineHS(config);
    CL.Preprocess(); 
    
%     %**********************************************
%     % Equilibration from off-equilibrium IC
%     %**********************************************
    rhoGas = CL.optsPhys.rhoGas_sat;
    rhoLiq = CL.optsPhys.rhoLiq_sat;

    [om,rho1D_wl,params] = CL.Compute1D('WL');            
	[om,rho1D_wg,params] = CL.Compute1D('WG');
    [om,rho1D_lg,params] = CL.Compute1D('LG');
    
    rho1D_wl = repmat(rho1D_wl,CL.IDC.N1,1);
    rho1D_wg = repmat(rho1D_wg,CL.IDC.N1,1);
    rho1D_lg = kronecker(rho1D_lg,ones(CL.IDC.N2,1));
    
    rho_ic = rho1D_wg + (rho1D_wl - rho1D_wg).*(rho1D_lg-rhoGas)/(rhoLiq - rhoGas);
    CL.x_eq = CL.optsPhys.kBT*log(rho_ic) + CL.Vext;            
    CL.ComputeDynamics();
    CL.PostprocessDynamics([4,5.5]);
    
       
    CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});
    %CL.PlotDynamicValue({'UV_t','entropy'},{'save','MovingFrameOfReference'});
    %CL.ComputeEquilibrium(struct('solver','Newton'));    
    %CL.ComputeEquilibrium();              
    

%     %*************************************************************
%     % Increasing/Decreasing attraction of wall - Spontaneous spreading
%     %*************************************************************
%     CL.ComputeEquilibrium(struct('solver','Picard'));
% 
%      CL.optsPhys.V1.tau            = 5;    
%      CL.optsPhys.V1.epsilon_w_end  = 0.856;%epsilon_w_end
%      CL.optsNum.plotTimes.t_int(2) = 200;
% %     
% %     CL.optsPhys                  = rmfield(CL.optsPhys,'BCWall_U');
%      CL.ComputeDynamics();
%      CL.PostprocessDynamics();
%      CL.PlotDynamicValue({'UV_t','entropy'},{'save'});
%     
%     %***********************            
%     CL.optsPhys.V1.epsilon_w_end  = 0.5;%epsilon_w_end    
%     CL.ComputeDynamics();
%     CL.PostprocessDynamics();
%     CL.PlotDynamicValue({'UV_t','entropy'},{'save'});

%     %**********************************************
%     % Moving the wall - forced wetting
%     %**********************************************
%     
%     CL.optsNum.plotTimes.t_int(2) = 10;
%     CL.ComputeDynamics();            
%     CL.PostprocessDynamics();
%     CL.PlotDynamicValue({'UV_t','entropy'},{'save'});    
%     
%     %***********************        
%     CL.optsPhys.BCWall_U.u_max = -0.2;
%     CL.ComputeDynamics();
%     CL.PostprocessDynamics();
%     CL.PlotDynamicValue({'UV_t','entropy'},{'save'});

%     
%     %**********************************************
%     % Equilibration from off-equilibrium IC
%     %**********************************************
%     CL.optsNum.plotTimes.t_int(2) = 100;
%     CL.optsPhys.V1             = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.94);
%     CL.optsPhys.ModifyEq_to_IC = struct('Mode','expY2','a0',3,'a1',3);
%     CL.ComputeDynamics();
%     CL.PostprocessDynamics();
%     CL.PlotDynamicValue({'UV_t','entropy'},{'save'});
%     
%     %***********************
%     CL.optsPhys.ModifyEq_to_IC.a0 = -3;
%     CL.ComputeDynamics();
%     CL.PostprocessDynamics();    
%     CL.PlotDynamicValue({'UV_t','entropy'},{'save'});
                
%end