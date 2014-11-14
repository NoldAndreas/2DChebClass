function data = SimulationDisk(N1,N2,L1,L2,vext)

    disp('** Simulation Disk **');

    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;    
    
    %Initialization
    if(nargin == 0)
        R       = 3;
        vext    = @Vext2;
        N1      = 20;
        N2      = 20;
        N       = [N1;N2];
        Origin  = [3,-2];
    else
        R = L1;
    end
    
    DC                 = Disc(v2struct(R,N,Origin));
    [Pts,Diff,Int,Ind] = DC.ComputeAll();
    Interp = DC.ComputeInterpolationMatrix((0:0.02:1)',(0:0.02:1)',true,true);
	Conv   = DC.ComputeConvolutionMatrix(@f2);        
    
    intBound         = struct('y1_l',DC.PhysSpace1(0),...
                              'y1_u',DC.PhysSpace1(1),...
                              'y2_l',DC.PhysSpace2(0),...
                              'y2_u',DC.PhysSpace2(1));
                       
    [V,Vdiff,VInt] = vext(Pts.y1_kv,Pts.y2_kv,intBound,'polar');    
    VP             = vext(Interp.pts1,Interp.pts2);                  
                      
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrors(vplot,VP,V,Vdiff,Diff);
    
    %Check Interpolation
    data.InterPol = max(abs(Interp.InterPol*V - VP));       
    display([' Error in Interpolation: ', num2str(data.InterPol)]);

    DC.plot(V);    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);
        
    %Check Convolution
    fP_Conv  = Conv*f1(Pts.y1_kv,Pts.y2_kv);
%    data.Conv = max(abs(Interp.InterPol*fP_Conv - fConv(Interp.pts1,Interp.pts2)));
%    display([' Error in Convolution: ', num2str(data.Conv)]);        
    
    data.N1 = N1; data.N2 = N2;            
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    DC.plot(V); 
    title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);
    DC.plot(fP_Conv);     
    title('Convolution');
    pbaspect([1 1 1]);
    
    %DoPerformanceTest(DC);
    
    %***************************************************************
    %Check Sphere Integration
    %***************************************************************
    N            = [14,14];
    shape        = v2struct(R,N);
    shape.sphere = true;
    DC           = Disc(shape);
    DC.ComputeIntegrationVector();
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(r,t)
        x = r.*cos(t);
        y = r.*sin(t);
        
        x0 = 1;
        rn = sqrt((x-x0).^2+y.^2);
        
        s = 1;
        z = 1/(s*sqrt(2*pi))*exp(-(rn.^2)/(2*s^2));  
        z(r == inf) = 0;
    end
    function z = f2(r,t)       
        s = 2;
        z = 1/(s*sqrt(2*pi))*exp(-(r.^2)/(2*s^2));
        z(r == inf) = 0;        
    end

 
end