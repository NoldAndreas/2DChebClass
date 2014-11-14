function data = SimulationInfiniteCapillary(N1,N2,L1,L2,vext)

    disp('** Simulation Infinite Capillary **');

    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 = 20;  L1    = 2; 
        N2 = 20;  L2    = 1;
        N  = [N1,N2];
        
        Conv = struct('L1',L1,'L2',L2,'N',[N1,N2]);
        
        vext  = @Vext7;
        a  = -5; b = 5;
        y2Min = 0; y2Max = 1;
    end       
    a = - inf; b = inf;    
    
    
    IC     = InfCapillary(v2struct(L1,y2Min,y2Max,N,Conv));
    x1PAbs = IC.CompSpace1(3);
    [Pts,Diff,Int,Ind] = IC.ComputeAll();    
    Interp = IC.ComputeInterpolationMatrix((-x1PAbs:0.02:x1PAbs)',(-1:0.02:1)',true,true);
    Conv   = IC.ComputeConvolutionMatrix(@f2,true);                         
    
    intBound = struct('y1_l',IC.PhysSpace1(-1),...
                      'y1_u',IC.PhysSpace1(1),...
                      'y2_l',IC.PhysSpace2(-1),...
                      'y2_u',IC.PhysSpace2(1));
                       

    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound,'cart');    
    [VP]             = vext(Interp.pts1,Interp.pts2);                          
    
    %Check Differentiation
    vplot     = Interp.InterPol*V;        
    data      = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                

    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(2,1,1);
    IC.plot(V,'SC');    title('Interpolation');    
        
    subplot(2,1,2);
    IC.plot(Conv(:,Pts.y1_kv < inf)*V(Pts.y1_kv < inf),'SC');
    title('Convolution');    
                                    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f2(r)%x,y)
        z = exp(-r.^2);%(x.^2+y.^2));
    end
end