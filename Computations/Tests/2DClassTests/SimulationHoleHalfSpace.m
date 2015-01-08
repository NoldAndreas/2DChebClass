function [HHS,res] = SimulationHoleHalfSpace()

    disp('** Simulation Hole Half Space **');
    AddPaths();    
    close all;
    
    %Initialization    
    N1 =  20;   N2 = 20;
    R   = 1;
    L   = 4;        
    OriginHole  = [2;5];
    N           = [N1;N2];        
    y2Wall      = 0;
    
    HHS  = HoleInHalfSpace(v2struct(OriginHole,R,y2Wall,N,L));
    Int  = HHS.ComputeIntegrationVector(); 
    
    [V,VInt] = VTest(HHS.Pts.y1_kv,HHS.Pts.y2_kv);
    
    figure;
    HHS.PlotGrid(); hold on;
    HHS.PlotGridLines();
    xlim([-6 10]);
    ylim([0 10]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
    PrintErrorPos(Int*V-VInt,'Integration');   
    
    m = N1*N2;
    PrintErrorPos(Int(1:m)*V(1:m)-pi/4/R^2,'Integration Upper HS'); 
    PrintErrorPos(Int(1+m:2*m)*V(1+m:2*m)-pi/4/R^2,'Integration Left ');
    PrintErrorPos(Int(1+2*m:3*m)*V(1+2*m:3*m)-pi/4/R^2,'Integration Right ');
    PrintErrorPos(Int(1+3*m:end)*V(1+3*m:end)-(pi/4*(R^(-2) - (y2Wall-OriginHole(2))^(-2) )),'Integration Lower ');
    
    
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VInt] = VTest(y1,y2)       
            V     = ((y1-OriginHole(1)).^2+(y2-OriginHole(2)).^2).^(-2);  
                                                    
            VInt  = 3/4*pi/(R^2) +...
                    pi/4*(R^(-2) - (y2Wall-OriginHole(2))^(-2) );
            
    end

end