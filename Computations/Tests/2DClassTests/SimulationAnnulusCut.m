function [WDG,res] = SimulationAnnulusCut

    disp('** Simulation AnnulusCut **');    	
    close all;
        
    Geometry = struct('R_in',1,'R_out',5,...
                      'h',4.5,'N',[30,30]);
        
    WDG                = AnnulusCut(Geometry);    
    WDG.ComputeIntegrationVector();
    
    figure;
    WDG.PlotGridLines();  hold on;
    WDG.PlotGrid();
   % xlim([-5 5]);
	hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);    
    res.fig_handles{1} = gcf;    
    
end