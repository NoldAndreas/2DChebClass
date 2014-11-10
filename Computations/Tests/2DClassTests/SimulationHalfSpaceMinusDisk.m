function SimulationHalfSpaceMinusDisk()

    disp('** Simulation HalfSpace Minus Disk **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  15;   N2 = 15;
    R       = 1;
    L1      = 1;    
    Origin  = [0.;0.5];
    y2Wall  = 0;
    N       = [N1;N2];
    Rmax    = Inf;

    SMD                = HalfSpaceMinusDisk(v2struct(Origin,R,N,L1,Rmax,y2Wall));
    Int                = SMD.ComputeIntegrationVector();
    
    SMD.PlotGridLines();    
    SMD.PlotGrid();
    xlim([(Origin(1)-8) (Origin(1)+8)]);
    ylim([y2Wall (Origin(2) + 20)]);
    
    PtsCart = SMD.GetCartPts();
        
    [V,VInt]   = VTest(PtsCart.y1_kv,PtsCart.y2_kv);   
           
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
    xlim([(Origin(1)-5) (Origin(1)+5)]);
    ylim([y2Wall (Origin(2) + 5)]);
    pbaspect([10 (Origin(2) + 5 - y2Wall) 5]);
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************             
    function [V,VInt] = VTest(y1,y2)       
            V     = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(-2);  
            
            VInt  = 1/2*pi/(R^2);            
            if((Origin(2)-y2Wall)>R)
                VInt = VInt + 1/4*pi/(R^2) + pi/4*(R^(-2) - (y2Wall-Origin(2))^(-2) );
            else
                Y    = (Origin(2)-y2Wall);
                VInt = VInt + ...
                    2*(0.2e1 * atan(sqrt(R ^ 2 - Y ^ 2) / Y) * R ^ 4 + 0.4e1 * asin(Y / R) * R ^ 2 * Y ^ 2 + 0.2e1 * (R ^ 2 - Y ^ 2) ^ (0.3e1 / 0.2e1) * Y - pi * R ^ 4 + 0.2e1 * sqrt(R ^ 2 - Y ^ 2) * Y ^ 3) / Y ^ 2 / R ^ 4 / 0.8e1;
                
                %VInt  = 3/4*pi/(R^2) +...
%                    pi/4*(R^(-2) - (y2Wall-Origin(2))^(-2) );
            end                        
    end

end