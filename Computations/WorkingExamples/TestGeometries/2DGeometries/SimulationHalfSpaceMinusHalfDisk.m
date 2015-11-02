function [SMD,res] = SimulationHalfSpaceMinusHalfDisk()

    disp('** Simulation HalfSpace Minus HalfDisk **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;
    R       = 1;
    L1      = 1;    
    Origin  = [0;0];
    N       = [N1;N2];
    Rmax    = Inf;

    SMD                = HalfSpaceMinusHalfDisk(v2struct(Origin,R,N,L1,Rmax));
    [Pts,Diff,Int,Ind] = SMD.ComputeAll();    
    Interp             = SMD.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
    
    figure;
    SMD.PlotGridLines();  hold on;
    SMD.PlotGrid();
    SMD.PlotIsoline(0,'y1');
    xlim([(Origin(1)-2) (Origin(1)+2)]);
    ylim([(Origin(2)-0.5) (Origin(2) + 3)]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
    PtsCart = SMD.GetCartPts();
        
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
    SMD.plot(V,'SC'); 
    title('Interpolation');
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
   function [V,VInt] = VTest(y1,y2)       
            V     = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(-2);  
            VInt  = pi/2*(1/R^2 - 1/Rmax^2);
   end

end