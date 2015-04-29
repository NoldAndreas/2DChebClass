function [data,res] = SimulationHalfSpace_Composed()%N1,N2,L1,L2,vext)

    disp('** Simulation Half Space Composed**');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)
        N      = [20;20];
        N2bound = 10;
        L1 = 2;    L2 = 2;        
        y2Min = -0.5; h = 1;        
        vext  = @Vext15;       
        alpha = pi/2;%4;
    end                
    
    HSC         = HalfSpace_Composed(v2struct(L1,L2,N,N2bound,y2Min,h,alpha));            
    HSCCartPts  = HSC.GetCartPts();
    V           = vext(HSCCartPts.y1_kv,HSCCartPts.y2_kv);
    
    figure('color','white','Position',[0 0 600 600]);
    HSC.PlotGridLines(); 
    HSC.PlotGrid();
    
    HSC.Sub_HalfSpace.PlotIsoline(0,'y2');
    HSC.PlotIsoline(-1/sqrt(2),'y1');
    HSC.PlotIsoline(1/sqrt(2),'y1');
	set(gca,'fontsize',20); set(gca,'linewidth',1.5);
    ylim([-0.5 10]); xlim([-0.5 0.5]*10);        	
    xlabel('$y_1$','Interpreter','Latex','fontsize',25);
    ylabel('$y_2$','Interpreter','Latex','fontsize',25);
    %SaveFigure('HalfSpaceGrid_Composed');
    res.fig_handles{1} = gcf;
    
    figure;
    HSC.plot(V,'SC');%,'SC'); 
    
    PhysAreaBX  = struct('y1Min',-3*L1,'y1Max',3*L1,...
                        'y2Min',0,'y2Max',10,...
                        'N',[40,40]);
    BX          = Box(PhysAreaBX);
    BXCartPts   = BX.GetCartPts();    
    
    
    Interp      = HSC.SubShapePtsCart(BX.Pts);    
    VP          = vext(BXCartPts.y1_kv,BXCartPts.y2_kv);

    vplot       = Interp*V;    
    data        = displayErrorsPos(HSC.Pts,vplot,VP,V);
    
end