function [data,res] = SimulationSlitMinusDisk(N1,N2,vext)

    disp('** Simulation Slit Minus Disk **');
    AddPaths();    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  20;   N2 = 20;        
        R       = 3;
        L1      = 2;
        y2Min   = -2.5;         
        y2Max   = 2; 
        OriginDisk  = [0;0];
        N       = [N1;N2];
        LeftRight = 'Left';        
    end        
    
    
    SMD                = SlitMinusDisk(v2struct(OriginDisk,R,y2Min,y2Max,N,L1,LeftRight));
    [Pts,Diff,Int,Ind] = SMD.ComputeAll();    
    Interp             = SMD.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
    
    figure;
    SMD.PlotGridLines();  hold on;
    SMD.PlotGrid();
    xlim([-5 5]);
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);    
    res.fig_handles{1} = gcf;    
        
    [V,Vdiff,VInt]   = VTest(Pts.y1_kv,Pts.y2_kv);
    [VP]             = VTest(Interp.pts1,Interp.pts2);
           
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
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
    function [V,VDiff,VInt] = VTest(y1,y2)       
            V           = y1.^(-4);  
            dVy1        = -4*y1.^(-5);
            dVy2        = zeros(size(y1));  
            %ddVy1       = -y2.*cos(y1);
            %ddVy2       = zeros(size(y1));
            %dVy1dVy2    = -sin(y1);

            VDiff       = struct('dy1',dVy1,'dy2',dVy2);
                                 %'ddy1',ddVy1,'ddy2',ddVy2,...
                                 %'dy1dy2',dVy1dVy2);            
                                                    
            VInt        = 1/(3*R^2)*( y2Max.*(R^2-y2Max.^2).^(-1/2) ...
                                    - y2Min.*(R^2-y2Min.^2).^(-1/2) );                                     
            
    end

end