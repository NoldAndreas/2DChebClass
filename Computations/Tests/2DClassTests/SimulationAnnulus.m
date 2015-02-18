function [data,res] = SimulationAnnulus()

    disp('** Simulation Annulus **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;        
    RMin     = 1;
    RMax     = 5;
    Origin  = [-2;3];
    N       = [N1;N2];    

    AS                 = Annulus(v2struct(Origin,RMin,RMax,N));
    [Pts,Diff,Int,Ind] = AS.ComputeAll();    
    Interp             = AS.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    
    
    figure;
    AS.PlotGridLines();    
    AS.PlotGrid();
    	
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
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
            V     = BarkerHendersonCutoff_2D(r);
            VInt  = pi*(1/RMin^2 - 1/RMax^2);
   end

end