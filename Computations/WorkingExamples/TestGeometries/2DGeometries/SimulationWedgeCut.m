function [WDG,res] = SimulationWedgeCut

    disp('** Simulation WedgeCut **');    
	if(length(dbstack) == 1)
        AddPaths();
    end        
    close all;
        
    Geometry = struct('R_in',1,'R_out',2,...
                      'h',1.5,'N',[10,10]);
    
    vext              = @Vext7;
   
    tic
    WDG                = WedgeCut(Geometry);
    [Pts,Diff,Int,Ind] = WDG.ComputeAll();
    Interp = WDG.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);    
    toc
    
    figure;
    WDG.PlotGridLines();  hold on;
    WDG.PlotGrid();
    xlim([-1.5 1.5]);
    ylim([-1.6 -0.7]);
	hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);    
    res.fig_handles{1} = gcf;    
    
     %Check Polar Spectral/Spectral map in 2D:                                             
     [V,Vdiff] = vext(Pts.y1_kv,Pts.y2_kv);   
     VP        = vext(Interp.pts1,Interp.pts2);           
       
%     figure;     
%     WDG.plot(V);
	 
     vplot     = Interp.InterPol*V;
     displayErrors(vplot,VP,V,Vdiff,Diff);
    
end