function ContactLineDynamics_45_90_135(kBT)
    
    AddPaths();

    if(nargin < 1)
        %kBT = 0.75;
        kBT = 0.9;
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
                   
	plotTimes = struct('t_int',[0,250],'t_n',100);
    
    PlotAreaCart = struct('y1Min',-7.5,'y1Max',7.5,...
                          'y2Min',0.5,'y2Max',15.5,...
                          'N1',100,'N2',100,'NFlux',10);    

    optsNum = struct('PhysArea',PhysArea,... 
                     'PlotAreaCart',PlotAreaCart,...                     
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
                           'zetaLiq',1.3,'zetaVap',0.01);
    %previous version
    %optsViscosity = struct('etaLiq',5,'etaVap',1,...
    %                           'zetaLiq',5,'zetaVap',1);
    %******************************************************
    %******************************************************
                       
    %BCwall        = struct('bc','sinHalf','tau',1);
	%BCwall        = struct('bc','exp','tau',1,'u_max',0.2);
    BCwall = [];

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',kBT,...                                               
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

    alphas  = [0,30,40,45,50,60,70,80,90,100,110,120,130,135,140,150];    
    eps     = DataStorage([],@FindEps,struct('alpha_degs',alphas,'config',config),[],[],{});
    
    %***************************************************
    optsDefault = {'snapshots'};
    dynRes = {};
    dynRes{end+1} = MCL_Comp(config,90,eps(alphas==60),optsDefault);    
    dynRes{end+1} = MCL_Comp(config,90,eps(alphas==120),optsDefault);
                                            
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==90),optsDefault);
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==0),optsDefault);                                            
                                            
    dynRes{end+1} = MCL_Comp(config,120,eps(alphas==90),optsDefault);
    
    %****************************************************
    % further computations
    optsDefault = {};
	dynRes{end+1} = MCL_Comp(config,90,eps(alphas==70),optsDefault);
    dynRes{end+1} = MCL_Comp(config,90,eps(alphas==80),optsDefault);
	dynRes{end+1} = MCL_Comp(config,90,eps(alphas==100),optsDefault);
    dynRes{end+1} = MCL_Comp(config,90,eps(alphas==110),optsDefault);
    
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==40),optsDefault);
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==50),optsDefault);
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==70),optsDefault);
    dynRes{end+1} = MCL_Comp(config,60,eps(alphas==80),optsDefault);
    
    dynRes{end+1} = MCL_Comp(config,120,eps(alphas==100),optsDefault);
    dynRes{end+1} = MCL_Comp(config,120,eps(alphas==110),optsDefault);
    dynRes{end+1} = MCL_Comp(config,120,eps(alphas==130),optsDefault);
    dynRes{end+1} = MCL_Comp(config,120,eps(alphas==140),optsDefault);
    
    dynRes = PostProcess(dynRes);
    GetDataForThesis(dynRes);
    
	PlotVelocitiesForThesis(dynRes,'contactlineVel_y1_0','$U_{\text{CL}}$');
    PlotVelocitiesForThesis(dynRes,'contactangle_0','$\theta$'); 
	PlotMKTComparison(dynRes);
    
    %***************************************************
     
%     MCL_Comp(config,{'90','advancing','snapshots'});
%     %MCL_Comp(config,{'45','advancing','snapshots'});        
%     %MCL_Comp(config,{'135','advancing','snapshots'});
%     
%     config.optsNum.plotTimes.t_int(2) = 200;
%     %MCL_Comp(config,{'90','receding','snapshots'});   
%     MCL_Comp(config,{'45','receding','snapshots'});        
%     MCL_Comp(config,{'135','receding','snapshots'});

    function eps = FindEps(input,misc)        
        input.config.optsNum.PhysArea.alpha_deg = 90;        
        eps  = FindEpwFromContactAngle(input.config,input.alpha_degs);
    end    
    function res = MCL_Comp(config,alpha_deg,epw,opts)
        
        recomp = true;%[];
         
        config.optsNum.PhysArea.alpha_deg = alpha_deg;
        config.optsPhys.V1.epsilon_w      = epw;
         
        if(alpha_deg <= 80)
            y1_Int = [-5,10];
        elseif(alpha_deg >= 100)
            y1_Int = [-10,5];
        else
            y1_Int = [-7.5,7.5];
        end

        config.optsNum.PlotAreaCart.y1Min = y1_Int(1);
        config.optsNum.PlotAreaCart.y1Max = y1_Int(2);
         
        input.config = config;
        input.opts   = opts;
        if(config.optsPhys.kBT == 0.75)
            input.fitInterval = [4 5.5];
        elseif(config.optsPhys.kBT == 0.9)
            input.fitInterval = [5.5 7];
        end
        
        res = DataStorage('MovingContactAngleResults',...
                        @ContactLineDynamics_X_degrees,...
                          input,{},recomp,{'fitInterval'}); %Eq: 170 degrees      
    end
    function GetDataForThesis(resAll)
             
        i_plot = 100;
        for i = 1:length(resAll)
            res = resAll{i};
            disp(['Data for ',...
                  num2str(res.thetaInitial),' -> ',...
                  num2str(res.thetaEq),...     
                  ' (theta_in -> theta_eq) ',...
                  ' for epw = ',num2str(res.epw),...
                  ' at t = ',num2str(res.t(i_plot)),...
                  ' : U_CL = ',num2str(res.contactlineVel_y1_0(i_plot)),...
                  ' theta = ',num2str(res.contactangle_0(i_plot))]);                
        end
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
        zeta_receding  = 4;%wRecSum/URecSum;
        zeta_advancing = 10;%wAdvSum/UAdvSum;
%         zeta_receding  = zeta_receding  / n_receding;
%         zeta_advancing = zeta_advancing / n_advancing;
        
        plot([0;0.07],[0;0.07]/zeta_receding,'k');
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
            plot(res{ii}.t, res{ii}.(varName),[res{ii}.col]); hold on; %res{ii}.lin
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
                            
            res{i}.col = cols{1+mod(i,length(cols))};
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
    
end


























































































































































































































































































































































































 















































































