function [fContour] =  PlotContourResults(this,options) %plain
    if(nargin == 1) %select all
        options = {'hI','hII','hIII','hI_alignedwith_hII','5contours',...
                            'YoungCA','hContour','Isoline','save','newFigure'};        
    elseif((nargin >= 2) && ~(iscell(options)))
        options = {options};
    end
    
    deltaY2 = -0.5;
    
    if(IsOption(options,'newFigure'))
        fContour = figure('Color','white','Position',[0 0 1200 800]);
    end
    
    rho           = GetRhoEq(this);
    rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
	rhoGas_sat    = this.optsPhys.rhoGas_sat;
    R             = this.optsPhys.sigmaS/2;            
    if(~isempty(this.y1_SpectralLine))
        y1  = this.y1_SpectralLine.Pts.y;
    end
            
	%*******************************************************
    % ******************** Contour Plot ********************
    %*******************************************************    
    if(IsOption(options,'Isoline'))
        if(~isempty(this.hContour))            
            plot(y1,this.hContour,'b-.'); %,'linewidth',2.5
        end
    end
    if(IsOption(options,'Isoline'))
        if(~isempty(this.IsolineInterfaceY2))
            plot(this.IsolineInterfaceY2,this.y2,'m-.','linewidth',1.5);
        end
    end
    if(IsOption(options,'YoungCA'))
        %Plot Young Contact Angle     
        y1MaxP = max(this.IDC.Interp.pts1);
        y2MaxP = max(this.IDC.Interp.pts2);

        theta_CS     = this.optsNum.PhysArea.alpha_deg*pi/180;
        alphaYoungIn = this.alpha_YCA;

        pt2.y1_kv = min(y1MaxP,y2MaxP/tan(theta_CS));
        pt2.y2_kv = tan(theta_CS)*pt2.y1_kv;

        dx1       = pt2.y1_kv - pt2.y2_kv/tan(alphaYoungIn);
        x1M       = min(y1MaxP - dx1,y2MaxP/tan(alphaYoungIn));
        %x1M = min(PlotArea.y1Max,PlotArea.y2Max/tan(alphaYoungIn));
        plot([dx1;x1M+dx1],[0;x1M*tan(alphaYoungIn)],'--m','linewidth',2); hold on;

        x1M = min(y1MaxP,y2MaxP/tan(alphaYoungIn));
        plot([0,x1M],[0,x1M*tan(theta_CS)],'--g','linewidth',2); hold on;
    end
    
    if(IsOption(options,'5contours'))
        optDetails.clabel = true;        
        %optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
        optDetails.nContours = [0.1,0.3,0.5,0.6,0.7];        
        this.IDC.plot(rho,'contour',optDetails);  hold on;                               	
    else
        optDetails.y2CartShift = deltaY2;
        optDetails.clabel      = false;  
        optDetails.linewidth   = 1.0;  
        
        PlotDensityContours(this,rho,optDetails);                         
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
    if(IsOption(options,'hIII') && ~isempty(this.hIII))% && ((nargin == 1) || ~plain))
        plot(y1,this.hIII+R+deltaY2,'k'); %'linewidth',1.5 adsorption/(rhoLiq_sat-rhoGas_sat)%%%%
        h0 = min(this.hIII);
    else
        h0 = 0;
    end

    if(IsOption(options,'hI_alignedwith_hII') && ~isempty(this.hI))        
        [DeltaY1_II,DeltaY1_III] = this.ComputeDeltaFit();                
        plot(this.y1_I+DeltaY1_II,this.hI+R+h0+deltaY2,'k-.'); %,'linewidth',1.5
    end
    
	if(IsOption(options,'hI') && ~isempty(this.hI))        
        [DeltaY1_II,DeltaY1_III] = this.ComputeDeltaFit();                
        plot(this.y1_I+DeltaY1_III,this.hI+R+h0+deltaY2,'k-.'); %,'linewidth',1.5
    end
    
    if(IsOption(options,'hII') && ~isempty(this.hII))
        plot(y1,this.hII+R+h0+deltaY2,'k--'); %'linewidth',1.5
    end  
    
    if(IsOption(options,'hIV') && ~isempty(this.hIV))
        plot(y1,this.hIV+R+h0+deltaY2,'r--');  %,'linewidth',1.5
    end  
    
    
	xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
	ylabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    
 %   pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 1]);
    if(IsOption(options,'save'))  
        SaveCurrentFigure(this,'contour');        
    end   
  
end