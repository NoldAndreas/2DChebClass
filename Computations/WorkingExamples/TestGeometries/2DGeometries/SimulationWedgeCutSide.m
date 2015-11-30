function [WDG,res] = SimulationWedgeCutSide

    disp('** Simulation WedgeCutSide **');    
	if(length(dbstack) == 1)
        AddPaths();
    end        
    close all;
        
    Geometry = struct('R_in',1,'R_out',2.5,...
                      'h',0.9,'N',[10,10],...
                      'leftRight','left');
    
    vext              = @Vext7;
   
    tic
    WDG                = WedgeCutSide(Geometry);
    [Pts,Diff,Int,Ind] = WDG.ComputeAll();
    Interp = WDG.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);    
    toc
    
    figure('color','white','Position',[0 0 600 300]);
    WDG.PlotGridLines();  hold on;
    WDG.PlotGrid();
    WDG.PlotIsoline(0,'y1');
   % xlim([-5 5]);
	hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);    
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);
        
    ylim([-1 -0.3]); xlim([-2.5 -0.3]);
    ax = get(gcf,'children');
    xl = get(ax,'xlim');
    yl = get(ax,'ylim');
	pbaspect(ax,[(xl(2)-xl(1)) (yl(2)-yl(1)) 1]);
    
    %pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1])
    %SaveFigure('WedgeCutSide');
    res.fig_handles{1} = gcf;    
    
     %Check Polar Spectral/Spectral map in 2D:                                             
     [V,Vdiff] = vext(Pts.y1_kv,Pts.y2_kv);   
     VP        = vext(Interp.pts1,Interp.pts2);           
       
%      figure;     
%      WDG.plot(V);
	 
     vplot     = Interp.InterPol*V;
     displayErrors(vplot,VP,V,Vdiff,Diff);
    
end