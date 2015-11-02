function [SMD,res] = SimulationHalfSpaceMinusDisk()

    disp('** Simulation HalfSpace Minus Disk **');
    AddPaths();    
    close all;
    
    %Initialization
    N1 =  20;   N2 = 20;
    R       = 1;
    L1      = 0.97;
    Origin  = [0.;7];%0.5]
    y2Wall  = 0.5;
    N       = [N1;N2];    

    SMD                = HalfSpaceMinusDisk(v2struct(Origin,R,N,L1,y2Wall));
    Int                = SMD.ComputeIntegrationVector();
    
    figure;
    SMD.PlotGridLines();    
    SMD.PlotGrid();
    xlim([(Origin(1)-3) (Origin(1)+3)]);
    ylim([y2Wall (Origin(2) + 4)]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
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
    SMD.plot(V,'SC'); 
    title('Interpolation');
    xlim([(Origin(1)-5) (Origin(1)+5)]);
    ylim([y2Wall (Origin(2) + 5)]);
    pbaspect([10 (Origin(2) + 5 - y2Wall) 1]);
    
                                                    
    %******* Check Interpolation *******
    for i = 1:2        
        InterpCPta      = SMD.SubShape{i}.GetCartPts(...
                                        SMD.SubShape{i}.Interp.pts1,...
                                        SMD.SubShape{i}.Interp.pts2);
        
        Vi           = VTest2(SMD.SubShape{i}.GetCartPts);
        VInterp      = SMD.SubShape{i}.Interp.InterPol*Vi;
        VInterpExact = VTest2(InterpCPta);
        PrintErrorPos(VInterp-VInterpExact,['Interpolation of subspace ',num2str(i)],InterpCPta);
    end
    [V,VInt] = VTest2(SMD.GetCartPts);
    PrintErrorPos(Int*V-VInt,'Integration of Barker Henderson');
    
 
    
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


    function [V,VInt] = VTest2(pts)
        y1 = pts.y1_kv;
        y2 = pts.y2_kv;
        d  = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(1/2);        
        V  = BarkerHenderson_2D(d);
        
        z    = Origin(2)-y2Wall;
        VInt = 2/3*pi/z^3 - 4/45*pi/z^9 - 25/32*pi^2;
    end

end