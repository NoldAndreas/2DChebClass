function [HS,res] = SimulationFittedHalfSpace()
    vext  = @Vext15;        
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,...
                       'L1',3,'L2',2,'N',[20;20],...
                       'y2Min',0);

    PlotArea       = struct('y1Min',-5,'y1Max',5,'L1',3,'L2',2,...
                       'N2',100,'N1',100,'y2Min',0,'y2Max',4);    

    Phys_Area.Conv      = struct('L',3,'L2',1,'N',[20,20]);%'ep2Conv',0.1     
    
    HS                             = HalfSpace(Phys_Area);%v2struct(L1,L2,N));
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(PlotArea);        
    
    figure;
    HS.PlotSubGrid(0,inf,@Phi2DLongRange);
    xlim([-5 5]);
    ylim([-5 5]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
    figure;
    HS.PlotSubGrid(0,10,@Phi2DLongRange);
    xlim([-5 5]);
    ylim([-10 5]);
	hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{2} = gcf;
    
    figure;
    HS.PlotSubGrid(0,0.5,@Phi2DLongRange);
    xlim([-4 4]);
    ylim([-1 7]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{3} = gcf;    

end