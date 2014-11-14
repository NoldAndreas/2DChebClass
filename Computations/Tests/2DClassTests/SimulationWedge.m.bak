function SimulationWedge

    disp('** SimulationPolarWedge **');    
	if(length(dbstack) == 1)
        AddPaths();
    end        
    close all;
    
    %Initialization
    
    R                   = 2;
    half_wedge_angle    = pi/2;%4*pi/5;
    Geometry = struct('R_in',0,'R_out',R,...
                      'th1',-half_wedge_angle,...
                      'th2',half_wedge_angle,'N',[30,20]);
    
    vext              = @Vext7;
    

    tic
    WDG                = Wedge(Geometry);
    [Pts,Diff,Int,Ind] = WDG.ComputeAll();
    Interp = WDG.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    Conv   = WDG.ComputeConvolutionMatrix(@f2,true);            
    toc
    
     %Check Polar Spectral/Spectral map in 2D:                      
     intBound         = struct('y1_l',WDG.PhysSpace1(-1),...
                               'y1_u',WDG.PhysSpace1(1),...
                               'y2_l',WDG.PhysSpace2(-1),...
                               'y2_u',WDG.PhysSpace2(1));
                       
     [V,Vdiff,VInt] = vext(Pts.y1_kv,Pts.y2_kv,intBound,'polar');    
     VP             = vext(Interp.pts1,Interp.pts2);           
       
     %Check Interpolation
     WDG.doPlots(V);
    
	 %Check Differentiation
     vplot     = Interp.InterPol*V;
     displayErrors(vplot,VP,V,Vdiff,Diff);
    
	 %Check Integration
     display([' Error in Integration: ', num2str(Int*V-VInt)]);
    

     %Check Convolution
     figure
     fP_Conv  = Conv*f1(Pts.y1_kv,Pts.y2_kv);
     WDG.doPlots(fP_Conv);
                          
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