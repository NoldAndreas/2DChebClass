function ThesisNanoscale_Fig11_ComputeDynamicContactLines()

    AddPaths('ThesisNanoscale');            
    close all;
    
%     config = ThesisNanoscale_GetStandardConfig(90,[]);
%     epw = FindEpwFromContactAngle(config,180);
%     disp(epw);     
      res{1} = DataStorage('MovingContactAngleResults',...
                           @DoDynamicComputation,...
                           struct('alpha_deg',90,'epw',0.594,'maxT',400),{}); %Eq: 120 degrees
                     
      res{2} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',90,'epw',1.071,'maxT',400),{}); %Eq: 60 degrees                     
                     
      res{3} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',45,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                          
                     
      res{4} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                                               
      
      res{5} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',1.27,'maxT',400),{}); %Eq: 0 degrees                                                               
      
      res{6} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.0,'maxT',400),{}); %Eq: 180 degrees         
      res{6}.erratic = true;                     
                                              
      res{7} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',135,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                                               
                     
      res{8} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',120,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                                                                    

	
      
    col = {'sk','vb','sm','vg','sr','vy'};
	figure('color','white');
    index = 90;
    for i = 1:6
        if(isfield(res{i},'erratic'))
            continue;
        end
        %plot(cos(res{i}.thetaInitial*pi/180)-cos(res{i}.thetaEq),res{i}.contactangle_0(end),['s',col{i}]);  hold on;
        plot(cos(res{i}.thetaInitial*pi/180)-cos(res{i}.thetaEq),res{i}.contactlineVel_y1_0(index),[col{i}]);  hold on;
        %plot(cos(res{i}.thetaInitial*pi/180)-cos(res{i}.thetaEq),res{i}.contactangle_0(end)*180/pi,'sk');  hold on;
    end
    
    for i = 1:6
        if(isfield(res{i},'erratic'))
            continue;
        end
        res{i}.cosDifference = cos(res{i}.contactangle_0*pi/180)-cos(res{i}.thetaEq);
        subplot(3,1,1);
        plot(res{i}.t,res{i}.contactlineVel_y1_0,[col{i}]); hold on;
        subplot(3,1,2);
        plot(res{i}.t,res{i}.contactlinePos_y1_0,[col{i}]); hold on;
        subplot(3,1,3);
        plot(res{i}.cosDifference,res{i}.contactlineVel_y1_0,[col{i}]); hold on;
    end
	
    %Analyse Contact Line motion
	
    
%     DoDynamicComputation(90,1.154,400); %Eq: 45 degrees
%     DoDynamicComputation(90,0.453,400); %Eq: 135 degrees
%     
%     DoDynamicComputation(45,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(60,1.154,400); %Eq: 45 degrees
%     
%     DoDynamicComputation(120,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(135,0.453,400); %Eq: 135 degrees
          
    function dynRes = DoDynamicComputation(input,inputMisc)
        %elements of 'input': alpha_deg,epw,maxT        
        alpha_deg = input.alpha_deg;
        epw       = input.epw;
        maxT      = input.maxT;
        
        try 
            config = ThesisNanoscale_GetStandardConfig(alpha_deg,epw,maxT);
            config.optsNum.PlotAreaCart       = struct('y1Min',-7.5,'y1Max',7.5,...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'N1',100,'N2',100,'NFlux',20);

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
            
            %******
            dynRes.epw                 = CL.optsPhys.V1.epsilon_w;
            dynRes.thetaEq             = CL.alpha_YCA;
            dynRes.thetaInitial        = CL.optsNum.PhysArea.alpha_deg;
            dynRes.t                   = CL.dynamicsResult.t;
            dynRes.contactlinePos_y1_0 = CL.dynamicsResult.contactlinePos_y1_0;
            dynRes.contactlineVel_y1_0 = CL.dynamicsResult.contactlineVel_y1_0;
            dynRes.contactangle_0      = CL.dynamicsResult.contactangle_0.val;
            
%             %For Thesis %,'fittedInterface'
             CL.PlotDynamicValue({'rho_t','contactangle_0'},...
                 {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize'});
%             
             CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},...
                 {'save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize'}); %'save'
%             
%             %For Movie:
%             CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});
            close all;
        catch err
            disp('ERROR')
            rethrow(err);        
        end
    end
    

end