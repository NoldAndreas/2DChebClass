function ThesisNanoscale_Fig11_ComputeDynamicContactLines()

    AddPaths('ThesisNanoscale');            
    close all;
    
%     config = ThesisNanoscale_GetStandardConfig(90,[]);
%     epw = FindEpwFromContactAngle(config,180);
%     disp(epw);     
      recomp = false;      
      
    %********************************************   
    %*** 90 degrees initial contact angle
    %********************************************      
     res90{1} = DataStorage('MovingContactAngleResults',...
                           @DoDynamicComputation,...
                           struct('alpha_deg',90,'epw',0.594,'maxT',400),{},recomp); %Eq: 120 degrees                              
                     
                       
     res90{2} = DataStorage('MovingContactAngleResults',...
                           @DoDynamicComputation,...
                           struct('alpha_deg',90,'epw',0.7,'maxT',400),{}); %Eq: 108.3 degrees                       
                       
                       
     res90{3} = DataStorage('MovingContactAngleResults',...
                           @DoDynamicComputation,...
                           struct('alpha_deg',90,'epw',0.8,'maxT',400),{}); 
                       
     res90{4} = DataStorage('MovingContactAngleResults',...
                           @DoDynamicComputation,...
                           struct('alpha_deg',90,'epw',0.9,'maxT',400),{});                          
                     
     res90{5} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',90,'epw',1.071,'maxT',400),{},recomp); %Eq: 60 degrees                     
                     	                        
    %********************************************   
    %*** 60 degrees initial contact angle
    %********************************************
    res60{1} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',0.856,'maxT',400),{},recomp); %Eq: 90 degrees                                                               

    res60{2} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',0.9,'maxT',400),{});                                                                                    
    res60{2}.erratic = true;
                     
                     
    res60{3} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',1.0,'maxT',400),{});                                                                                       
    res60{3}.erratic = true;

    res60{4} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',1.1,'maxT',400),{});                                                                                       
                     
    res60{5} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',1.2,'maxT',400),{});                       
     
    res60{6} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',60,'epw',1.27,'maxT',400),{},recomp); %Eq: 0 degrees                                                               
                                                 
    %********************************************             
    %*** 60 degrees initial contact angle
    %********************************************             
    
    res120{1} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.1,'maxT',400),{},recomp); %Eq: 170 degrees        
    res120{1}.erratic = true;
                     
    res120{2} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.2,'maxT',400),{}); 
    res120{2}.erratic = true;

    res120{3} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.4,'maxT',400),{}); 
                       
    res120{4} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.6,'maxT',400),{}); 
    
    res120{5} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                           struct('alpha_deg',120,'epw',0.8,'maxT',400),{}); 

	res120{6} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',120,'epw',0.856,'maxT',400),{},recomp); %Eq: 90 degrees                                                                  
    
     %********************************************                  
    
      res{1} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',45,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                          
                                
                                              
      res{2} = DataStorage('MovingContactAngleResults',...
                         @DoDynamicComputation,...
                         struct('alpha_deg',135,'epw',0.856,'maxT',400),{}); %Eq: 90 degrees                                                               
      res{2}.erratic = true;            
                     

                       
   %  res{10} = DataStorage('MovingContactAngleResults',...
%                           @DoDynamicComputation,...
%                           struct('alpha_deg',90,'epw',1.0,'maxT',400),{}); %Eq: 70.86 degrees                                              

    %********************************************                                      
    % Output results for thesis
    %********************************************                                      
    
    %ewf = 1.071, theta_n = 90
    GetDataForThesis(res90{5});        
    %ewf = 0.594, theta_n = 90
    GetDataForThesis(res90{1});    
    %ewf = 1.270, theta_n = 60
    GetDataForThesis(res60{6});        
    %ewf = 0.856, theta_n = 60
    GetDataForThesis(res60{1});
    %ewf = 0.100, theta_n = 120
    GetDataForThesis(res120{1});    
    %ewf = 0.856, theta_n = 120
    GetDataForThesis(res120{6});    
    
    %********************************************                                      
    %********************************************                                      

    
    res60 = PostProcess(res60);    
    res120 = PostProcess(res120);    
    res90 = PostProcess(res90);    
    res = PostProcess(res);    
    
    res = [res90,res60,res120,res];    
    
    PlotMKTComparison(res);
	PlotVelocitiesForThesis({res90{1},res90{5}},'contactlineVel_y1_0','$U_{\text{CL}}$');
    PlotVelocitiesForThesis({res90{1},res90{5}},'contactangle_0','$\theta$');    
    
    
    f2 = figure('color','white','Position',[0 0 250 200]);
    for i = 6:2:8%16:length(res)
        if(isfield(res{i},'erratic'))
            continue;
        end        
        mark = (20:100);
        
        subplot(2,2,1);
        plot(res{i}.t, res{i}.contactlineVel_y1_0,[res{i}.lin,res{i}.col]); hold on;
        plot(res{i}.t, res{i}.contactlineVel_y1_0,['o',res{i}.lin,res{i}.col]); hold on;
        subplot(2,2,2);
        plot(res{i}.t,res{i}.contactangle_0,res{i}.col); hold on;
        plot(res{i}.t,res{i}.contactangle_0,['o',res{i}.col]); hold on;
        subplot(2,2,3);
        plot(res{i}.cosDiff_t(mark),res{i}.contactlineVel_y1_0(mark),res{i}.col); hold on;
        plot(res{i}.cosDiff_t(mark),res{i}.contactlineVel_y1_0(mark),['o',res{i}.col]); hold on;
        subplot(2,2,4);
        plot(res{i}.GDiff_t(mark),res{i}.contactlineVel_y1_0(mark),['o',res{i}.col]); hold on;
        plot(res{i}.GDiff_t(mark),res{i}.contactlineVel_y1_0(mark),res{i}.col); hold on;
    end
    pbaspect([1 1 1]);
    %xlabel('$t$','Interpreter','Latex');
    %ylabel('$U_{CL}$','Interpreter','Latex');
    
    %inset2(f1,f2,0.4,[0.3,0.5]);
    %close(f2);    
    SaveFigure('ContactLineMeasurement_All_APSDFD2015');
	
    %Analyse Contact Line motion
	
%     DoDynamicComputation(90,1.154,400); %Eq: 45 degrees
%     DoDynamicComputation(90,0.453,400); %Eq: 135 degrees
%     
%     DoDynamicComputation(45,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(60,1.154,400); %Eq: 45 degrees
%     
%     DoDynamicComputation(120,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(135,0.453,400); %Eq: 135 degrees

    function GetDataForThesis(res)
        i_plot = 100;
        disp(['Data for epw = ',num2str(res.epw),...
              ' theta_in = ',num2str(res.thetaInitial),...
              ' at t = ',num2str(res.t(i_plot)),...
              ' : U_CL = ',num2str(res.contactlineVel_y1_0(i_plot)),...
              ' theta = ',num2str(res.contactangle_0(i_plot))]);                
    end

    function PlotMKTComparison(res)
        gammaLV =  0.2853;
        
        figure('color','white','Position',[0 0 300 250]);
        index = 80;
        leg = {};
        
        %zeta_advancing = 0; n_advancing = 0;
        %zeta_receding  = 0; n_receding  = 0;        
        wAdvSum = 0; wRecSum = 0;
        UAdvSum = 0; URecSum = 0;
        
        for i = 1:length(res)
            if(isfield(res{i},'erratic'))
                continue;
            end

            %subplot(1,2,1);
            w = gammaLV*res{i}.cosDiff_t(index);
            U = res{i}.contactlineVel_y1_0(index);
            plot(w,U,res{i}.sym,'MarkerFaceColor',res{i}.col,'MarkerEdgeColor',res{i}.col);  hold on;
            %plot(res{i}.cosDiff,res{i}.contactlineVel_y1_0(index),...
            %subplot(1,2,2);
            %plot(res{i}.GDiff,res{i}.contactlineVel_y1_0(index),...
             %           res{i}.sym,'MarkerFaceColor',res{i}.col,'MarkerEdgeColor',res{i}.col);  hold on;
             
            leg = [leg,res{i}.text];
            
            if(U ~= 0)
                if(res{i}.thetaEq < res{i}.thetaInitial)
                    wAdvSum = wAdvSum + w;
                    UAdvSum = UAdvSum + U;
%                     n_advancing    = n_advancing + 1;
%                     zeta_advancing = zeta_advancing + w/U;
%                     disp(w/U);
                else
                    wRecSum = wRecSum + w;
                    URecSum = URecSum + U;
                    %n_receding     = n_receding + 1;
                    %zeta_receding  = zeta_receding + w/U;
                end
            end
            
        end
        zeta_receding  = 8;%wRecSum/URecSum;
        zeta_advancing = 30;%wAdvSum/UAdvSum;
%         zeta_receding  = zeta_receding  / n_receding;
%         zeta_advancing = zeta_advancing / n_advancing;
        
        plot([0;0.11],[0;0.13]/zeta_receding,'k');
        plot([0;-0.13],[0;-0.13]/zeta_advancing,'--k');
        
        pbaspect([1 1 1]);

        %subplot(1,2,1);
        xlabel('$\surfaceTensionLV\klamm{\cos(\theta)-\cos(\theta_{eq})}$','Interpreter','Latex');
        ylabel('$U_{\text{CL}}$','Interpreter','Latex');
        legend(leg,'Location','SouthOutside','Interpreter','Latex');
        %subplot(1,2,2);
        %xlabel('$G(\theta_{in})-G(\theta_{eq})$','Interpreter','Latex');
        %ylabel('$U_{CL}$','Interpreter','Latex');
        SaveFigure('ContactLineMeasurement_MKTComparison');
        disp(['zeta_advancing = ',num2str(zeta_advancing)]);
        disp(['zeta_receding = ',num2str(zeta_receding)]);
    end

    function PlotVelocitiesForThesis(res,varName,yLabel)
        figure('color','white','Position',[0 0 200 150]);
        for ii = 1:length(res)
            if(isfield(res{ii},'erratic'))
                continue;
            end        
            %mark = (20:100);            
            plot(res{ii}.t, res{ii}.(varName),[res{ii}.lin,res{ii}.col]); hold on;
           % plot(res{i}.t, res{i}.contactlineVel_y1_0,['o',res{i}.lin,res{i}.col]); hold on;            
        end
        pbaspect([1 1 1]);
        xlabel('$t$','Interpreter','Latex');
        ylabel(yLabel,'Interpreter','Latex');

        %inset2(f1,f2,0.4,[0.3,0.5]);
        %close(f2);    
        SaveFigure(['ContactLine_t_',varName]);
    end
          
    function res = PostProcess(res)
        
        cols = {'r','m','k','b','c','g'};
        
        for i = 1:length(res)        
            if(isfield(res{i},'erratic'))
                continue;
            end

            res{i}.thetaEq   = res{i}.thetaEq*180/pi;
            res{i}.cosDiff   = cos(res{i}.thetaInitial*pi/180)-cos(res{i}.thetaEq*pi/180);
            res{i}.cosDiff_t = cos(res{i}.contactangle_0*pi/180)-cos(res{i}.thetaEq*pi/180);
            %res{i}.cosDiff_t = cos(res{i}.thetaInitial*pi/180)-cos(res{i}.contactangle_0*pi/180);

            lambdaEta      = 0.2;
            res{i}.GDiff   = GHR_lambdaEta(res{i}.thetaInitial*pi/180,lambdaEta)-...
                             GHR_lambdaEta(res{i}.thetaEq*pi/180,lambdaEta);

            %res{i}.GDiff_t   = GHR_lambdaEta(res{i}.contactangle_0*pi/180,lambdaEta)-...
    %                           GHR_lambdaEta(res{i}.thetaEq*pi/180,lambdaEta);            

            res{i}.GDiff_t   = GHR_lambdaEta(res{i}.thetaInitial*pi/180,lambdaEta) - ...
                                GHR_lambdaEta(res{i}.contactangle_0*pi/180,lambdaEta);                                   
                            
            res{i}.col = cols{i};
%            if(res{i}.thetaEq < res{i}.thetaInitial)
%                 res.advancing = 
%                 res{i}.sym = 's';
%                 res{i}.lin = '-';
%             else
%                 res{i}.sym = 'v';
%                 res{i}.lin = '--';
%             end
            
            switch res{i}.thetaInitial
                case 45
                    res{i}.sym = 'v';
                case 60
                    res{i}.sym = 's';
                case 90 
                    res{i}.sym = 'o';
                case 120
                    res{i}.sym = '^';
                case 135
                    res{i}.sym = 'd';
            end
            
            res{i}.text = ['$',num2str(res{i}.thetaInitial),'^\circ \rightarrow ',num2str(res{i}.thetaEq,3),'^\circ$'];
        end

    end
    function dynRes = DoDynamicComputation(input,inputMisc)
        %elements of 'input': alpha_deg,epw,maxT        
        alpha_deg = input.alpha_deg;
        epw       = input.epw;
        maxT      = input.maxT;
        
        try 
            config = ThesisNanoscale_GetStandardConfig(alpha_deg,epw,maxT);
            config.optsNum.PlotAreaCart       = struct('y1Min',-7.5,'y1Max',7.5,...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'N1',100,'N2',100,'NFlux',10);

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
                 {'save','MovingFrameOfReference','Snapshots','dimensionlessLabels','start','contour','PublicationSize','CLy2Shift'});
%             
             CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},... %
                 {'save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
             
             %without entropy plot
             CL.PlotDynamicValue({'rho_t','fittedInterface','UV_t','contactangle_0'},... %
                 {'contour','save','MovingFrameOfReference','Snapshots','streamlines','dimensionlessLabels','end','PublicationSize','CLy2Shift'}); %'save'
%             
             disp(['final CL velocity: ',num2str(CL.dynamicsResult.contactlineVel_y1_0(end))]);
%             %For Movie:
%             CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});
            close all;
        catch err
            disp('ERROR')
            rethrow(err);        
        end
    end
    

end