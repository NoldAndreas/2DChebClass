function dynRes = ContactLineDynamics_X_degrees(input,misc) %config,alpha_deg,epw,opts)

    config    = input.config;    
    opts      = input.opts;    

    AddPaths('CodePaper');            
    close all;        
    
    %************************************************
    %************************************************

    CL = ContactLineHS(config);
    CL.Preprocess(); 

    %**********************************************
    % Equilibration from off-equilibrium IC
    %**********************************************

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
    CL.PostprocessDynamics([4 5.5]);
    CL.PlotInterfaceFittingQuality([25,50,75,100]);
    %CL.PostprocessDynamics();

    %********************************************************************
    dynRes.epw                 = CL.optsPhys.V1.epsilon_w;
    dynRes.thetaEq             = CL.alpha_YCA;
    dynRes.thetaInitial        = CL.optsNum.PhysArea.alpha_deg;
    dynRes.t                   = CL.dynamicsResult.t;
    dynRes.contactlinePos_y1_0 = CL.dynamicsResult.contactlinePos_y1_0;
    dynRes.contactlineVel_y1_0 = CL.dynamicsResult.contactlineVel_y1_0;
    dynRes.contactangle_0      = CL.dynamicsResult.contactangle_0.val;
    %********************************************************************
    
    if(IsOption(opts,'videos'))
        CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});        
    end
    
    if(IsOption(opts,'snapshots'))  
             CL.PlotDynamicValue({'rho_t','contactangle_0','pathlines'},...
                         {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize','CLy2Shift'});
                     
             CL.PlotDynamicValue({'rho_t','contactangle_0','streaklines'},...
                         {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize','CLy2Shift'});                     
                     
             CL.PlotDynamicValue({'rho_t','contactangle_0','pathlines_MovingFrameOfReference'},...
                         {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize','CLy2Shift'});                     
        
             CL.PlotDynamicValue({'rho_t','contactangle_0','pathlines_MovingFrameOfReference'},...
                 {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize','CLy2Shift'});
             
             CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},... %
                 {'save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
                          
             CL.PlotDynamicValue({'rho_t','fittedInterface','UV_t','contactangle_0',},... %
                 {'contour','save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
             
             disp(['final CL velocity: ',num2str(CL.dynamicsResult.contactlineVel_y1_0(end))]);

            CL.optsNum.PlotAreaCart       = struct('y1Min',-5,'y1Max',5,...
                                                   'y2Min',0.5,'y2Max',7.5,...
                                                   'N1',40,'N2',40,'NFlux',10);
            CL.IDC.InterpolationPlotCart(CL.optsNum.PlotAreaCart,true);
            figure('color','white','Position',[0 0 300 300]);                        
            CL.PlotDynamicValue({'entropy','fittedInterface'},... 
                 {'3D','save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'

            figure('color','white','Position',[0 0 300 300]);                        
            CL.PlotDynamicValue({'comprEntr','fittedInterface'},... 
                 {'3D','save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
             
            figure('color','white','Position',[0 0 300 300]);                        
            CL.PlotDynamicValue({'shearEntr','fittedInterface'},... 
                 {'3D','save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
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
%
%
%******************************************************************
%******************************************************************
%******************************************************************
%******************************************************************
%******************************************************************
%******************************************************************
%******************************************************************
%******************************************************************
%
%     if(nargin < 1)
%         opts = {'advancing'};
%     end   
%     
%     if(IsOption(opts,'135'))
%         if(IsOption(opts,'advancing'))
%             config.optsNum.PhysArea.alpha_deg = 135;
%             config.optsPhys.V1.epsilon_w       = 0.856; %= 90 degree contact angle
%         elseif(IsOption(opts,'receding'))
%             config.optsNum.PhysArea.alpha_deg = 120;
%             config.optsPhys.V1.epsilon_w       = 0.453; %= 135 degree contact angle
%         else
%             return;
%         end
%         config.optsNum.PlotAreaCart =     struct('y1Min',-10,'y1Max',5,...
%                                                  'y2Min',0.5,'y2Max',15.5,...
%                                                  'N1',100,'N2',100,'NFlux',10);                               
%     elseif(IsOption(opts,'90'))        
%         config.optsNum.PhysArea.alpha_deg = 90;
%         if(config.optsPhys.kBT == 0.75)
%             if(IsOption(opts,'advancing'))                
%                 %epw       = 1.3; % = +/- 0 degree contact angle
%                 config.optsPhys.V1.epsilon_w  =  1.154; %= 45 degree contact angle
%             elseif(IsOption(opts,'receding'))    
%                 %epw       = 0.05; %= 180 degree contact angle
%                 config.optsPhys.V1.epsilon_w  = 0.453; %= 135 degree contact angle        
%             end    
%         else
%             if(IsOption(opts,'advancing'))
%                 config.optsPhys.V1.epsilon_w = FindEpwFromContactAngle(config,45);
%             elseif(IsOption(opts,'receding'))
%                 config.optsPhys.V1.epsilon_w = FindEpwFromContactAngle(config,135);        
%             else
%                 config.optsPhys.V1.epsilon_w = FindEpwFromContactAngle(config,90);
%             end
%         end                
%         config.optsNum.PlotAreaCart = struct('y1Min',-7.5,'y1Max',7.5,...
%                                              'y2Min',0.5,'y2Max',15.5,...
%                                              'N1',100,'N2',100,'NFlux',10);
%     elseif(IsOption(opts,'45'))
%         if(IsOption(opts,'receding'))
%             config.optsNum.PhysArea.alpha_deg = 45;
%             config.optsPhys.V1.epsilon_w      = 0.856; %= 90 degree contact angle        
%         elseif(IsOption(opts,'advancing'))
%             config.optsNum.PhysArea.alpha_deg = 60;
%             config.optsPhys.V1.epsilon_w      = 1.22;  %= +/- 30 degree contact angle
%             %epw       = 1.154; %= 45 degree contact angle        
%         else
%             return;
%         end
%         config.optsNum.PlotAreaCart = struct('y1Min',-5,'y1Max',10,...
%                                              'y2Min',0.5,'y2Max',15.5,...
%                                              'N1',100,'N2',100,'NFlux',10);
%         config.optsNum.V2Num.N       = [80,80];        
%         config.optsNum.FexNum.N1disc = 50;
%         config.optsNum.FexNum.N2disc = 50;  
%     end     
end