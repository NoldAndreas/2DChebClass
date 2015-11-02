function [data,res] = SimulationBox_M1(N1,N2,L1,L2,vext)

    disp('** Simulation Box M1 **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  25;   N2 = 20; N = [N1;N2];
        
        y1Min = -1; y1Max = 1;
        y2Min = -1; y2Max = 1;
                
        vext  = @Vext12;
    end        
    
    BI     = Box_Interface(v2struct(y1Min,y1Max,y2Min,y2Max,N));
    [Pts,Diff,Int,Ind] = BI.ComputeAll();
    Interp             = BI.ComputeInterpolationMatrix((-1:0.02:1)',...
                                (-1:0.02:1)',true,true);            
            
    intBound = struct('y1_l',-1,'y1_u',1,'y2_l',-1,'y2_u',1); 
 
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound);    
    [VP]             = vext(Interp.pts1,Interp.pts2);                          
    
    BI.plot(V,'SC');    
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    subplot(2,2,1); BI.plot(Vdiff.ddy2,'SC');
    subplot(2,2,2); BI.plot(Diff.DDy2*V,'SC');
    subplot(2,2,3); BI.plot(Diff.DDy2*V-Vdiff.ddy2,'SC');
    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
        
    %Check Convolution
    %fP = f1(Pts.y1_kv,Pts.y2_kv);
    %data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    %display([' Error in Convolution: ', num2str(data.Conv)]);
    
    data.N1 = N1; data.N2 = N2;
%    data.Interp = Interp; data.f = V;    
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    BI.plot(V);    
    title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);
    BI.plot(fConv(Interp.pts1,Interp.pts2));
    title('Convolution');
    pbaspect([1 1 1]);
    
    figure;
    BI.PlotGridLines();    
    BI.PlotGrid();
    	
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    %***************************************************************
    %   Mapping functions:
    %***************************************************************         
    function [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_PhysLinear(x1,x2)
        %[y1_kv,y2_kv,J,DDy1,DDy2] = Comp_to_Phys(x1,x2)
        %[y1_kv,dy1dx1] = LinearMap(x1,-L1,L1); %,dx1,ddx1,dddx,ddddx]                
        n  = length(x1);

        [y1_kv,dy1dx1] =  LinearMap(x1,-L1,L1);
        [y2_kv,dy2dx2] =  LinearMap(x2,-L2,L2);
        
        if(nargout >= 3)
            J        = zeros(n,2,2);
            J(:,1,1) = dy1dx1;
            J(:,2,2) = dy2dx2;
        end
        
        if(nargout >= 4)
            dH1        = zeros(n,2,2);
        end
        
        if(nargout >= 4)
            dH2        = zeros(n,2,2);
        end
        
        %[z,dz,dx,ddx,dddx,ddddx] = QuadMapAB(xR,L1,-L1,L1);
    end
%     
%     
%     function [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_Phys(x1,x2)
%         %[y1_kv,y2_kv,J,DDy1,DDy2] = Comp_to_Phys(x1,x2)
%         %[y1_kv,dy1dx1] = LinearMap(x1,-L1,L1); %,dx1,ddx1,dddx,ddddx]                
%         n  = length(x1);
%         c1 = 0.4;
%         c3 = 0.05;
% 
%         [y1_kv,Diffy1] =  M1Tref(x1,x2*c1+c3,0.07);                 
%         [y2_kv,dy2dx2] =  LinearMap(x2,-L2,L2);
%         
%         if(nargout >= 3)
%             J        = zeros(n,2,2);
%             J(:,1,1) = Diffy1.dydx; 
%             J(:,1,2) = Diffy1.dyd_d*c1;           
%             J(:,2,1) = zeros(n,1);
%             J(:,2,2) = dy2dx2;            
%         end
%         
%         if(nargout >= 4)
%             dH1      = zeros(n,2,2);
%             dH1(:,1,1) = Diffy1.dyddx; %d2y1dx1dx1; 
%             dH1(:,1,2) = c1*Diffy1.dydxd_d;%d2y1dx1dx2;             
%             dH1(:,2,1) = c1*Diffy1.dydxd_d;%d2y1dx1dx2;     
%             dH1(:,2,2) = (c1^2).*Diffy1.dydd_d;
%         end
%         
%         if(nargout >= 4)
%             dH2        = zeros(n,2,2);
%             dH2(:,1,1) = zeros(n,1); %d2y2dx1dx1; 
%             dH2(:,1,2) = zeros(n,1); %d2y2dx1dx2;             
%             dH2(:,2,1) = zeros(n,1); %d2y2dx1dx2;     
%             dH2(:,2,2) = zeros(n,1); %d2y2dx2dx2;             
%         end
%         
%         %[z,dz,dx,ddx,dddx,ddddx] = QuadMapAB(xR,L1,-L1,L1);
%     end
%     function [x1,x2] = Phys_to_Comp(y1,y2)
%          %xf  = InvQuadMapAB(z,L1,-L1,L1);
%          x1  = y1/L1;   x2  = y2/L2;
%     end    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(x,y)
        z = x;
    end
    function z = f2(x,y)
        z = exp(-(x.^2+y.^2));
    end
    function z = fConv(x,y)
        z = (-sqrt(pi)*exp(-(x+1).^2) - pi*x.*erf(1+x) +...
              sqrt(pi)*exp(-(x-1).^2) + pi*x.*erf(-1+x)).*...
              (-erf(1+y)+erf(y-1))/4;
    end
end