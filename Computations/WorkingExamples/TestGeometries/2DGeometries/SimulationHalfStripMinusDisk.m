function [SMD,res] = SimulationHalfStripMinusDisk()

    disp('** Simulation Strip Minus Disk **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;        
    R       = 1;
    L1      = 1;
    y2Wall  = -0.9;    
    Origin  = [0;0];
    N       = [N1;N2];
    LeftRight = 'Right';%'Left';
    TopBottom = 'Bottom';    

    SMD                = HalfStripMinusDisk(v2struct(Origin,R,y2Wall,N,L1,LeftRight,TopBottom));
    [Pts,Diff,Int,Ind] = SMD.ComputeAll();    
    Interp             = SMD.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
    
    figure;%('color','white','Position',[0 0 600 600]);
    SMD.PlotGridLines();  hold on;
    SMD.PlotGrid();
    SMD.PlotIsoline(0,'y1');
    if(strcmp(LeftRight,'Left'))
        xlim([-5 0]);
    else
        xlim([0 3]);
    end
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
            VInt  = 1/4*pi/(R^2) +...
                    pi/4*(R^(-2) - (y2Wall-Origin(2))^(-2) );             
   end

end