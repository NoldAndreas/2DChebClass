function SimulationHalfStripMinusDisk()

    disp('** Simulation Strip Minus Disk **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;        
    R       = 3;
    L1      = 0.5;
    y2Wall  = 0;    
    Origin  = [0;3];
    N       = [N1;N2];
    LeftRight = 'Left';%'Left';
    %Rmax      = 5;


    SMD                = HalfStripMinusDisk(v2struct(Origin,R,y2Wall,N,L1,LeftRight));
    [Pts,Diff,Int,Ind] = SMD.ComputeAll();    
    Interp             = SMD.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
    
    SMD.PlotGridLines();  hold on;
    SMD.PlotGrid();
    if(strcmp(LeftRight,'Left'))
        xlim([-5 0]);
    else
        xlim([0 5]);
    end
    
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
    SMD.doPlots(V,'SC'); 
    title('Interpolation');
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
   function [V,VInt] = VTest(y1,y2)       
            V     = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(-2);  
            VInt  = 1/4*pi/(R^2) +...
                    pi/4*(R^(-2) - (y2Wall-Origin(2))^(-2) );             
   end

end