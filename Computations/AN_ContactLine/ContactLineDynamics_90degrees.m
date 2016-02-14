function ContactLineDynamics_90degrees(config,opts)

    if(nargin < 1)
        opts = {'advancing'}; %'chemical'
    end    

    AddPaths('CodePaper');            
    close all;
    
    %advancing = false;
    alpha_deg = 90;
   
%     if(IsOption(opts,'dragging'))        
%             epw      = 0.856; %= 90 degree contact angle
%             maxT     = 20;            
%         if(IsOption(opts,'advancing'))        
%             BCwall   = struct('bc','exp','tau',1,'u_max',-0.2);            
%         elseif(IsOption(opts,'receding'))                
%             BCwall   = struct('bc','exp','tau',1,'u_max',0.2);            
%         end
%     elseif(IsOption(opts,'chemical'))    
    if(IsOption(opts,'advancing'))                
        %epw       = 1.3; % = +/- 0 degree contact angle
        epw       =  1.154; %= 45 degree contact angle
    elseif(IsOption(opts,'receding'))    
        %epw       = 0.05; %= 180 degree contact angle
        epw       = 0.453; %= 135 degree contact angle        
    end
%    end
    
	PlotAreaCart = struct('y1Min',-7.5,'y1Max',7.5,...
                          'y2Min',0.5,'y2Max',15.5,...
                          'N1',100,'N2',100,'NFlux',20);
 
    config.optsNum.PlotAreaCart       = PlotAreaCart;
    config.optsNum.PhysArea.alpha_deg = alpha_deg;
    config.optsPhys.V1.epsilon_w      = epw;
        
    CL = ContactLineHS(config);
    CL.Preprocess(); 
    
    if(IsOption(opts,'dragging'))
        CL.ComputeEquilibrium(struct('solver','Picard'));
    else    
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
    end
    CL.ComputeDynamics();
    CL.PostprocessDynamics([4 5.5]);
    CL.PlotInterfaceFittingQuality([1,26,51,76,100]);
    
    if(IsOption(opts,'videos'))
        CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});
    end
    %CL.PlotDynamicValue({'UV_t','entropy'},{'save'});
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