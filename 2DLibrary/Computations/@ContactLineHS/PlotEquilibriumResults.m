function [fContour] =  PlotEquilibriumResults(this,plain,saveFigs)
    %Compute contact angle from density profile
  %  [alphaM,pt1,pt2] = MeasureContactAngle();    
    
    %Compute excess grand potential as a function of y1
    %Fex_Y1 = GetExcessGrandPotential_Y1(rho);        
    
    %**************************************
    %Initialization
    
    global dirData    
    
%     if((nargin>2) && ~isempty(bounds1))                        
%         PlotArea.y1Min = bounds1(1);
%         PlotArea.y1Max = bounds1(2);
%         PlotArea.y2Min = bounds2(1);
%         PlotArea.y2Max = bounds2(2);
%         
%         PlotArea.N1    = 100;
%         PlotArea.N2    = 100;
%         
%         InitInterpolation(this,PlotArea);
%     elseif(isempty(this.IDC.Interp))
%         InitInterpolation(this);
%     end
    
    %PlotArea = this.optsNum.PlotArea;
    
    %PlotArea.y1Min = min(this.IDC.Interp.pts1);
    %PlotArea.y1Max = max(this.IDC.Interp.pts1);
    %PlotArea.y2Min = min(this.IDC.Interp.pts2);
    %PlotArea.y2Max = max(this.IDC.Interp.pts2);
    
    rho           = GetRhoEq(this);
    rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
	rhoGas_sat    = this.optsPhys.rhoGas_sat;
    R             = this.optsPhys.sigmaS/2;        
    y1            = this.y1_SpectralLine.Pts.y;
            
	%*******************************************************
    % ******************** Contour Plot ********************
    %*******************************************************
    if((nargin < 5) || (saveFigs))
        fContour = figure('Color','white','Position',[0 0 1200 800]);
    end
    
    if((nargin < 5) || ~plain)
        optDetails.clabel = true;        
        %optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
        optDetails.nContours = [0.1,0.3,0.5,0.6,0.7];        
        this.IDC.doPlots(rho,'contour',optDetails);  hold on;  
        
        %Plot Young Contact Angle     
        y1MaxP = max(this.IDC.Interp.pts1);
        y2MaxP = max(this.IDC.Interp.pts2);

        theta_CS  = this.optsNum.PhysArea.alpha_deg*pi/180;
        alphaYoungIn = this.alpha_YCA;

        pt2.y1_kv = min(y1MaxP,y2MaxP/tan(theta_CS));
        pt2.y2_kv = tan(theta_CS)*pt2.y1_kv;

        dx1       = pt2.y1_kv - pt2.y2_kv/tan(alphaYoungIn);
        x1M       = min(y1MaxP - dx1,y2MaxP/tan(alphaYoungIn));
        %x1M = min(PlotArea.y1Max,PlotArea.y2Max/tan(alphaYoungIn));
        plot([dx1;x1M+dx1],[0;x1M*tan(alphaYoungIn)],'--m','linewidth',2);            

        x1M = min(y1MaxP,y2MaxP/tan(alphaYoungIn));
        plot([0,x1M],[0,x1M*tan(theta_CS)],'--g','linewidth',2);
        
     	%**********************************************************
        % Add Isoline Thickness plot   
        if(~isempty(this.hContour))
            plot(y1,this.hContour,'b-.','linewidth',2.5);
        end
        
        if(~isempty(this.IsolineInterfaceY2))
            plot(this.IsolineInterfaceY2,this.y2,'m-.','linewidth',2.5);
        end
    else
        optDetails.y2CartShift = -0.5;
        optDetails.clabel = false;  
        optDetails.linewidth = 1.4;  
        
        %optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
        drho = rhoLiq_sat - rhoGas_sat;
        
        optDetails.nContours = rhoGas_sat + 0.1*drho;
        optDetails.linecolor = 'b';
        optDetails.linestyle = '--';
        this.IDC.doPlots(rho,'contour',optDetails);  hold on;  
        
        optDetails.nContours = rhoGas_sat + 0.5*drho;
        optDetails.linecolor = [0 0.75 0];
        this.IDC.doPlots(rho,'contour',optDetails);  hold on;  
        
        optDetails.nContours = rhoGas_sat + 0.9*drho;
        optDetails.linecolor = 'r';
        this.IDC.doPlots(rho,'contour',optDetails);  hold on;  
        
        xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
        ylabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    end    
    %Plot Measured Contact Angle
    %plot([pt1.y1_kv;pt2.y1_kv],[pt1.y2_kv;pt2.y2_kv],'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');        
    %plot((y2I-b)/slope,y2I,'--k','linewidth',1.5);   
    
    % Add Film Thickness plot   
%     shapeBox       = this.optsNum.PlotArea;
%     shapeBox.y2Max = this.optsNum.PlotArea.y2Max + 4;
%     shapeBox.N     = [40,120];
%     BX             = Box(shapeBox);
%     IP_BX          = this.IDC.SubShapePtsCart(BX.Pts);
%     [h1s,h2s,Int2] = BX.ComputeIntegrationVector();
%     Int2BX         = kronecker(eye(shapeBox.N(1)),Int2);
% 
%     %rho_wg_ref     = kronecker(ones(N1,1),rho1D_wg); 
%     adsorption     = Int2BX*IP_BX*(rho - rhoGas_sat);%rho_wg_ref);
%     hold on;     
   
    if(~isempty(this.hIII) && ((nargin == 1) || ~plain))
        plot(y1,this.hIII+R,'k-.','linewidth',2.5); %adsorption/(rhoLiq_sat-rhoGas_sat)
    end
    
    if(~isempty(this.hII) && ((nargin == 1) || ~plain))
        plot(y1,this.hII+R,'k--','linewidth',2.5);
    else

 %   pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 1]);
    if((nargin < 5) || saveFigs)      
        print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_contour'],gcf);
        saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_contour.fig']);
    end
    
    if((nargin >= 5) && plain)
        return;
    end
    
    %*******************************************************
    % ******************** 3D Plot ********************
    %*******************************************************
    figure('Color','white','Position',[0 0 1200 800]);
    this.IDC.doPlots(rho,'SC');
    zlabel('$\varrho$','Interpreter','Latex','fontsize',26);
    colormap(hsv);
    set(gca, 'CLim', [0, 1.0]);
    PlotArea = this.optsNum.PlotAreaCart;
    pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 5]);
    view([-10 5 3]);

    if((nargin < 5) || saveFigs)
        print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq],gcf);
        saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '.fig']);
    end
    %*******************************************************
    % ***************** Interface Plots ********************
    %*******************************************************
    figure('Color','white','Position',[0 0 800 1000],'name','1D Interface Plots');
    subplot(3,1,1);
    this.IDC.do1DPlotParallel(this.rho1D_lg); 
    title('Liquid-Gas Interface');
    ylabel('$\varrho$','Interpreter','Latex'); 
    xlabel('$y_1$','Interpreter','Latex');

    subplot(3,1,2);
    this.IDC.do1DPlotNormal(this.rho1D_wg);
    title('Wall-Gas Interface');
    ylabel('$\varrho$','Interpreter','Latex');  
    xlabel('$y_2$','Interpreter','Latex');

    subplot(3,1,3);
    this.IDC.do1DPlotNormal(this.rho1D_wl);
    title('Wall-Liquid Interface');
    ylabel('$\varrho$','Interpreter','Latex'); 
    xlabel('$y_2$','Interpreter','Latex');

    if((nargin == 1) || saveFigs)
        print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Interfaces'],gcf);
        saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Interfaces.fig']);
    end            
end