function Post_HFrom2DDisjoiningPressure(this,f1)
           
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;    
        
    h0            = min(abs(this.filmThickness));
    
    ratio_y_x =  (this.optsNum.PlotArea.y2Max- this.optsNum.PlotArea.y2Min)/...
                    (this.optsNum.PlotArea.y1Max- this.optsNum.PlotArea.y1Min);
                
    LP_1 = this.optsNum.PlotArea.y1Max - this.optsNum.PlotArea.y1Min;
    LP_2 = this.optsNum.PlotArea.y2Max - this.optsNum.PlotArea.y2Min;
    LP   = 20;
    
	if(nargin < 2)
        f1 = figure('Color','white','Position',[0 0 800 600]);
    end
    
	optDetails.clabel = false;  
	optDetails.linecolor = 'k';
    drho = rhoLiq_sat - rhoGas_sat;
    
    optDetails.y2CartShift = -0.5;
    optDetails.linewidth = 1;
	optDetails.nContours = rhoGas_sat + [0.05,0.5,0.95]*drho;
	%this.HS.doPlots(this.rho_eq,'contour',optDetails);  hold on;  
    PlotEquilibriumResults(this,[],[],true,false);
        
%    plot(this.y1,h0 + h_2D,'k--','linewidth',1.5); hold on;    
    plot(this.y1_I,h0 + this.hI,'k-.','linewidth',1.5);
 %   plot(this.y1,this.filmThickness,'k','linewidth',1.5);               
        
	xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
	ylabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    xlim([this.optsNum.PlotArea.y1Min this.optsNum.PlotArea.y1Max]);   
    %pbaspect([2 1 1]);

    D = this.DiffY1;
    D2 = this.DiffYY1;            
	
end