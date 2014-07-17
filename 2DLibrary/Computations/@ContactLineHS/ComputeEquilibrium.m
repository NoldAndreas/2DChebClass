function sol = ComputeEquilibrium(this,redo)    

    global dirData
    fprintf('Solving for equilibrium condition...\n');    
    
    %*****************************
    %Initialization
    %*****************************
	rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
	rhoGas_sat    = this.optsPhys.rhoGas_sat;
    kBT           = this.optsPhys.kBT;
    PtsCart     = this.HS.GetCartPts();  
     
    %*****************************
    p         = (this.rho1D_lg-rhoGas_sat)/(rhoLiq_sat-rhoGas_sat);    
    rho_ig    = kron(p,this.rho1D_wl) + kron(1-p,this.rho1D_wg);         
    x_ig      = kBT*log(rho_ig)+this.Vext;
    
    opts             = this.optsNum.PhysArea;
    opts.optsPhys    = this.optsPhys;    
    if(isfield(opts.optsPhys,'gamma'))
        opts.optsPhys    = rmfield(opts.optsPhys,'gamma');
    end
    if(isfield(opts.optsPhys,'Inertial'))
        opts.optsPhys    = rmfield(opts.optsPhys,'Inertial');
    end
    
    if(isfield(opts.optsPhys.V1,'epsilon_w_end'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_end');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_Amplitude'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_Amplitude');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_max'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_max');                
    end
    if(isfield(opts.optsPhys.V1,'tau'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'tau');
    end

    opts.maxComp_y2  = this.optsNum.maxComp_y2;
    opts.Comments    = this.configName;
    
    mark             = (PtsCart.y2_kv <= this.optsNum.maxComp_y2);     
            
    %misc = v2struct(mark,Vext,VAdd,Conv,IntMatrFex,x_ig);
    misc.mark = mark;       misc.Vext       = this.Vext; 
    misc.VAdd = this.VAdd;
    misc.x_ig = x_ig;
    misc.Conv = this.Conv;  misc.IntMatrFex = this.IntMatrFex; 
    misc.FexNum = this.optsNum.FexNum;
    
    if(nargin == 1)
        [sol,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                            @ComputeEquilibriumCondition,opts,misc); %true      
    else
        [sol,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                            @ComputeEquilibriumCondition,opts,misc,redo); %true      
    end

    this.rho_eq      = sol.rho;
    this.FilenameEq  = paramsEq.Filename;
	sol.Filename     = [dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename];
    
    close all;
    
 
%     %Compute contact angle from density profile
%     rhoV = (rhoLiq_sat+rhoGas_sat)/2; x2 = 0;
%     [alphaM,pt1,pt2] = MeasureContactAngle();    
%     
%     %Compute excess grand potential as a function of y1
%     %Fex_Y1 = GetExcessGrandPotential_Y1(rho);
%     sol.Filename = [dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename];
%     if(PersonalUserOutput)        
%         
%         figure('Color','white','Position',[0 0 1200 800]);
%         optDetails.clabel = true;        
%         optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
%         HS.doPlots(rho,'contour',optDetails);  hold on;  
%         plot([pt1.y1_kv;pt2.y1_kv],[pt1.y2_kv;pt2.y2_kv],'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');        
%         plot((y2I-b)/slope,y2I,'--k','linewidth',1.5);
%         %Plot Young Contact Angle                
%         dx1 = pt2.y1_kv - pt2.y2_kv/tan(alphaYoungIn);
%         x1M = min(PlotArea.y1Max-dx1,PlotArea.y2Max/tan(alphaYoungIn));
%         plot([dx1;x1M+dx1],[0;x1M*tan(alphaYoungIn)],'--m','linewidth',2);        
%         pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 1]);
%         
%         x1M = min(PlotArea.y1Max,PlotArea.y2Max/tan(alphaYoungIn));
%         plot([0,x1M],[0,x1M*tan(PhysArea.alpha_deg*pi/180)],'--g','linewidth',2);
%         
%         %**********************************************************
%         % Add Film Thickness plot        
%         shapeBox       = PlotArea;
%         shapeBox.y2Max = PlotArea.y2Max + 4;
%         shapeBox.N     = [40,120];
%         BX             = Box(shapeBox);
%         IP_BX          = HS.SubShapePtsCart(BX.Pts);
%         [h1s,h2s,Int2] = BX.ComputeIntegrationVector();
%         Int2BX         = kronecker(eye(shapeBox.N(1)),Int2);
%         
%         %rho_wg_ref     = kronecker(ones(N1,1),rho1D_wg); 
%         adsorption     = Int2BX*IP_BX*(rho - rhoGas_sat);%rho_wg_ref);
%         hold on; 
%         plot(BX.Pts.y1,adsorption/(rhoLiq_sat-rhoGas_sat)+R,'k-.','linewidth',2.5);
%         
%         if(saveFigs)       
%             print2eps([dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_contour'],gcf);
%             saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_contour.fig']);
%         end
%         
%         %*******************************************************
%         % ******************** Contour Plot ********************
%         %*******************************************************
%         figure('Color','white','Position',[0 0 1200 800]);
%         HS.doPlots(rho,'SC');
%         zlabel('$\varrho$','Interpreter','Latex','fontsize',26);
%         colormap(hsv);
%         set(gca, 'CLim', [0, 1.0]);
%         pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 5]);
%         view([-10 5 3]);
%         
%         if(saveFigs)       
%             print2eps([dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename],gcf);
%             saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '.fig']);
%         end
%         
%         %*******************************************************
%         % ***************** Interface Plots ********************
%         %*******************************************************
%         figure('Color','white','Position',[0 0 800 1000],'name','1D Interface Plots');
%         subplot(3,1,1);
%         HS.do1DPlotParallel(rho1D_lg); 
%         title('Liquid-Gas Interface');
%         ylabel('$\varrho$','Interpreter','Latex');
%         xlabel('$y_1$','Interpreter','Latex');
%         
%         subplot(3,1,2);
%         HS.do1DPlotNormal(rho1D_wg);
%         title('Wall-Gas Interface');
%         ylabel('$\varrho$','Interpreter','Latex');
%         xlabel('$y_2$','Interpreter','Latex');
%         
%         subplot(3,1,3);
%         HS.do1DPlotNormal(rho1D_wl);
%         title('Wall-Liquid Interface');
%         ylabel('$\varrho$','Interpreter','Latex');
%         xlabel('$y_2$','Interpreter','Latex');
%         
%         if(saveFigs)
%             print2eps([dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_Interfaces'],gcf);
%             saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_Interfaces.fig']);
%         end        
%     end
%     
%     Fex = GetExcessGrandPotential_Y1(rho);
%     if(~isfield(optsNum,'plotTimes_T'))
%         sol.Measured_StaticCA = alphaM;        
%         diary                
%         return;
%     end
end