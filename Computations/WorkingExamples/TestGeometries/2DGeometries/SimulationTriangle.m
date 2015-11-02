function [data,res] = SimulationTriangle(N1,N2,L1,L2,vext)

    disp('** Simulation Triangle **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  20;   N2 = 20;
        L1 = 4; L2 = 5;
        Y = [-1,0;L1+2,2;L1,L2];     
        vext  = @VTest;
    end            
    
    TR = Triangle(v2struct(Y),[N1,N2]);
    [Pts,Diff,Int] = TR.ComputeAll();    
    Interp = TR.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv);
    [VP]             = vext(Interp.pts1,Interp.pts2);
        
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
        
    subplot(2,2,1); TR.plot(Vdiff.dy1,'SC');
    subplot(2,2,2); TR.plot(Diff.Dy1*V,'SC');
    subplot(2,2,3); TR.plot(Diff.Dy1*V-Vdiff.dy1,'SC');    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
            
    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color            
    TR.plot(V);    
    title('Interpolation');    
    pbaspect([1 1 1]);  
    
    figure;
    TR.PlotGridLines();  hold on;
    TR.PlotGrid();
   % xlim([-5 5]);
	hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);    
    res.fig_handles{1} = gcf;   
    
    %***************************************************************
    %   Mapping functions:
    %***************************************************************             
    function [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_Phys(x1,x2)        
        n  = length(x1);

        [y1_kv,dy1dx1] =  LinearMap(x1,0,L1);
        [y2_kv_h,dy2dx2] =  LinearMap(x2,0,L2);
        
        y2_kv = y2_kv_h.*y1_kv/L1;
        
        if(nargout >= 3)
            J        = zeros(n,2,2);
            J(:,1,1) = dy1dx1;
            J(:,1,2) = zeros(n,1); 
            J(:,2,1) = y2_kv_h.*dy1dx1/L1;
            J(:,2,2) = y1_kv.*dy2dx2/L1;            
        end
        
        if(nargout >= 4)
            dH1        = zeros(n,2,2);            
        end
        
        if(nargout >= 4)
            dH2        = zeros(n,2,2);            
            dH2(:,1,2) = dy1dx1.*dy2dx2/L1;
            dH2(:,2,1) = dy1dx1.*dy2dx2/L1;            
        end
               
    end
    function [x1,x2] = Phys_to_Comp(y1,y2)         
         x1  = -1 + 2*y1/L1;
         x2  = -1 + 2*y2/(L2*y1);
         x2(y1 == 0 && y2 == 0) = 0;
    end    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest(y1,y2)       
            V           = (y1.^2).*cos(y2);  
            dVy1        = 2*y1.*cos(y2);        
            dVy2        = -(y1.^2).*sin(y2);  
            ddVy1       = 2*cos(y2);        
            ddVy2       = -(y1.^2).*cos(y2);  
            dVy1dVy2    = -2*y1.*sin(y2);  

            VDiff       = struct('dy1',dVy1,'dy2',dVy2,...
                                'ddy1',ddVy1,'ddy2',ddVy2,...
                                'dy1dy2',dVy1dVy2);
            
            VInt        = ((2-L2^2)*cos(L2)+2*L2*sin(L2)-2)*L1^3/(L2^3);           
    end


end