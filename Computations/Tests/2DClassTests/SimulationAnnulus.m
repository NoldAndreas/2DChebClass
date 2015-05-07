function [data,res] = SimulationAnnulus()

    disp('** Simulation Annulus **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;
    RMin     = 1;
    RMax     = 2.55;
    Origin  = [0,0];%[-2;3];
    N       = [N1;N2];    

    AS                 = Annulus(v2struct(Origin,RMin,RMax,N));
    [Pts,Diff,Int,Ind] = AS.ComputeAll();    
    Interp             = AS.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    
    
    figure('color','white','Position',[0 0 600 600]);
    AS.PlotGridLines();    
    AS.PlotGrid();    	
    AS.PlotIsoline(0,'y1');
    
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);            
    
    set(gca,'fontsize',20); set(gca,'linewidth',1.5);
    xlim([-3 3]); ylim([-3 3]);
    
    ax = get(gcf,'children');
    xl = get(ax,'xlim');
    yl = get(ax,'ylim');
	pbaspect(ax,[(xl(2)-xl(1)) (yl(2)-yl(1)) 1]);
    
    %pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1])
    SaveFigure('FiniteAnnulus');
    
    res.fig_handles{1} = gcf;    
 
    Origin = [2;3];
    AS.Origin = Origin;
    
    PtsCart = AS.GetCartPts();    
        
    [V,VInt]   = VTest(PtsCart.y1_kv,PtsCart.y2_kv);
    VP         = VTest(Interp.pts1,Interp.pts2);
           
    %Check Differentiation
    %vplot   = Interp.InterPol*V;        
    %data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
        
    %******** Plotting **********
    figure('Color','white')
    AS.plot(V,'SC'); 
    title('Interpolation');
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
   function [V,VInt] = VTest(y1,y2)       
            r     = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(0.5);  
            %V     = r.^(-4);  
            global r_cutoff
            r_cutoff = RMax;
            rc       = r_cutoff;
            V     = BarkerHendersonCutoff_2D(r);
            VInt  = -pi * (-(3080 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 4) - (5619 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1)) + 0.6000e4 * atan(sqrt((rc ^ 2 - 1))) * (rc ^ 10) - 0.7400e4 * sqrt((rc ^ 2 - 1)) * (rc ^ 6) - (4326 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 2) + 0.525e3 * sqrt((rc ^ 2 - 1)) + (6000 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 6) + 0.525e3 * sqrt((rc ^ 2 - 1)) * (rc ^ 2) + 0.350e3 * sqrt((rc ^ 2 - 1)) * (rc ^ 4)) / (rc ^ 10) / 0.3840e4;
           
            
%            VInt  = pi*(1/RMin^2 - 1/RMax^2);
   end

end