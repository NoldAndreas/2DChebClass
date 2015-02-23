function [y2,theta] = PlotInterfaceAnalysisY2(this,yInt)    

    %********************************************************************
    %Film thickness            
    if((nargin == 1) && ~isempty(this.IsolineInterfaceY2))
       iso = this.IsolineInterfaceY2;
       D   = this.DiffY2;
       y2  = this.y2;         
    else
        [y2,~,D] = InitAnalysisGridY(this,yInt,400);
        iso      = ComputeInterfaceContourY2(this,0.5,y2);
    end
        
    f1 = figure('color','white','Position',[0 0 800 800]);       
    
    theta          = atan(1./(D*iso))*180/pi;
    theta(theta<0) = theta(theta<0) + 180;
    plot(y2,theta,'k','linewidth',1.5); hold on;
        
    plot([min(y2) max(y2)],[this.optsNum.PhysArea.alpha_deg this.optsNum.PhysArea.alpha_deg],'k--','linewidth',1.5);    
    plot([this.optsNum.maxComp_y2 this.optsNum.maxComp_y2],[min(theta) max(theta)],'k--','linewidth',1.5);
    
    %Plot Measured limits and values:
    plot([min(y2) max(y2)],[this.alpha_YCA*180/pi this.alpha_YCA*180/pi],'k:','linewidth',1.5);
    plot([min(this.y2) max(this.y2)],[this.CA_deg_measured this.CA_deg_measured],'k-.','linewidth',1.5);
    %plot([min(this.y2) min(this.y2)],[]);
    %plot([max(this.y2) max(this.y2)],[]);
      
    xlabel('$y$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',25);
    
    xlim([min(y2) max(y2)]);
    ylim([min(theta) max(theta)]);
    set(gca,'fontsize',20);   
    
    SaveCurrentFigure(this,'InterfaceSlopeY2');

%    f2 = PlotEquilibriumResults(this,[-5 10],[0.5 (this.optsNum.maxComp_y2+3)],true,true);
%    plot([1 2],[3 4]);
%    inset2(f1,f2,0.45,[0.25,0.1]);
%    close(f2);            
%    
%    print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_InterfaceSlopeY2'],f1);
%    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_InterfaceSlopeY2.fig']);
    
end
