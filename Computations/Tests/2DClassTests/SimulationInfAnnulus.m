function SimulationInfAnnulus()

    disp('** Simulation InfAnnulus **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;        
    RMin    = 1;
    L       = 1;
    Origin = [2,3];
    N       = [N1;N2];    

    AS                = InfAnnulus(v2struct(Origin,RMin,L,N));
    [Pts,Diff,Int,Ind] = AS.ComputeAll();    
    Interp             = AS.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
    
    AS.PlotGridLines();  hold on;
    AS.PlotGrid();
    
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
            V     = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(-2);  
            VInt  = pi/RMin^2;
   end

end